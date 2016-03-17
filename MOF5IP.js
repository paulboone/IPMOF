require('./MOF5_Extended.js');
require('./MOF5.js');
// Interpenetration Functions --------------------------------------------------
var E1 = [22.156, 52.873, 30.213, 62.441]; // Epsilon for H, C, O, Zn
var S1 = [2.571, 3.431, 3.118, 2.462];     // Sigma for H, C, O, Zn
var A1 = ['H', 'C', 'O', 'Zn'];
var UC1 = [26, 26, 26];
var difAtomsInfo = [A1, E1, S1];

// Trilinear Interpolation
var trInterpolate = function(coor, numAtoms){
	var coor1 = []; var coor0 = []; var dif = [];

	for(var i = 0; i < coor.length; i++){
		coor0[i] = Math.floor(coor[i]);
		coor1[i] = Math.ceil(coor[i]);
		dif[i] = ( coor[i] - coor0[i] ) / ( coor1[i] - coor0[i] );
	};

	var i000 = findEmapIndex(coor0, grid.decimal, minCoor, maxCoor);
	var i100 = findEmapIndex([coor1[0],coor0[1],coor0[2]], grid.decimal, minCoor, maxCoor);
	var i001 = findEmapIndex([coor0[0],coor0[1],coor1[2]], grid.decimal, minCoor, maxCoor);
	var i101 = findEmapIndex([coor1[0],coor0[1],coor1[2]], grid.decimal, minCoor, maxCoor);
	var i010 = findEmapIndex([coor0[0],coor1[1],coor0[2]], grid.decimal, minCoor, maxCoor);
	var i110 = findEmapIndex([coor1[0],coor1[1],coor0[2]], grid.decimal, minCoor, maxCoor);
	var i011 = findEmapIndex([coor0[0],coor1[1],coor1[2]], grid.decimal, minCoor, maxCoor);
	var i111 = findEmapIndex(coor1, grid.decimal, minCoor, maxCoor);

	var V = []; var c00, c01, c10, c11, c0, c1, c;
	for(var i = 0; i < numAtoms; i++){
		c00 = energyMap[i000][i+3]*(1-dif[0]) + energyMap[i100][i+3]*dif[0];
		c01 = energyMap[i001][i+3]*(1-dif[0]) + energyMap[i101][i+3]*dif[0];
		c10 = energyMap[i010][i+3]*(1-dif[0]) + energyMap[i110][i+3]*dif[0];
		c11 = energyMap[i011][i+3]*(1-dif[0]) + energyMap[i111][i+3]*dif[0];

		c0 = c00*(1-dif[1]) + c10*dif[1];
		c1 = c01*(1-dif[1]) + c11*dif[1];

		c = c0*(1-dif[2]) + c1*dif[2];
		V[i] = c;
	};
	return V;
};

var calculateEnergyMap = function(atomsInfo, UC, extendedStructure, cutOff){
	var A = atomsInfo[0]; var E = atomsInfo[1]; var S = atomsInfo[2];
	var V = [];
	var Vtotal = [];
	var i = 0; var r, atomCoor, eps, sig;
	var minA = Math.round(- UC[0] / 2); var maxA = -minA;
	var minB = Math.round(- UC[1] / 2); var maxB = -minB;
	var minC = Math.round(- UC[2] / 2); var maxC = -minC;

	for(var x = minA; x <= maxA; x = x + grid.length){
	  for(var y = minB; y <= maxB; y = y + grid.length){
	    for(var z = minC; z <= maxC; z = z + grid.length){
	      energyMap[i] = [x,y,z];
	      for(var v = 0; v < E.length; v++){ Vtotal[v] = 0; };
	      for(var j = 0; j < extendedStructure.length; j++){
	        atomCoor = [extendedStructure[j][0], extendedStructure[j][1], extendedStructure[j][2]];
	        r = coorDist(energyMap[i],atomCoor);
	        if(r > cutOff){ continue; };
	        if(r === 0){ for(var v = 0; v < E.length; v++){Vtotal[v] = Infinity} }
	        else {
	          for(var a = 0; a < E.length; a++){
							if(extendedStructure[j][3] === A[a]){
		            for(var b = 0; b < E.length; b++){
		              eps = Math.sqrt(E[a]*E[b]);
		              sig = (S[a] + S[b]) / 2;
		              V[b] = 4*eps*(Math.pow((sig/r),12)-Math.pow((sig/r),6));
		              Vtotal[b] = Vtotal[b] + V[b];
		            };
							};
	          };
	        };
	      };
				for(var v = 0; v < E.length; v++){ energyMap[i].push(Vtotal[v]) };
	      i++;
	    };
	  };
	};
}; // Calculte energy Map Function

var findEmapIndex = function(coor, decimal, minCoor, maxCoor){
	var gridSize = Math.abs(minCoor) + Math.abs(maxCoor);
	for(var i = 0; i < 3; i++){
		coor[i] = Math.round(coor[i]*decimal)/decimal;
	};
	var eMapIndex = (coor[0]-minCoor)*Math.pow(gridSize+1,2) + (coor[1]-minCoor)*(gridSize+1) + coor[2]-minCoor;
	return eMapIndex;
};

var PBC = function(coor, cutOff){
	var pbcCoor = []; var translationAmount;
	for(var i = 0; i < coor.length; i++){
		if(Math.abs(coor[i]) > cutOff){
			translationAmount = Math.ceil(Math.abs(coor[i] / cutOff));
			pbcCoor[i] = coor[i] - Math.sign(coor[i])*(translationAmount)*cutOff;
		} else{
			pbcCoor[i] = coor[i];
		}
	};
	return pbcCoor;
};

var findAtomType = function(atomSymbol){
	var atomIndex;
	switch(atomSymbol){
		case 'H':
			atomIndex = 3;
			break;
		case 'C':
			atomIndex = 4;
			break;
		case 'O':
			atomIndex = 5;
			break;
		case 'Zn':
			atomIndex = 6;
			break;
	};
	return atomIndex;
};

var coorDist = function(coor1, coor2){
  var dx = coor2[0] - coor1[0];
  var dy = coor2[1] - coor1[1];
  var dz = coor2[2] - coor1[2];
  var distance = Math.sqrt(Math.pow(dx,2) + Math.pow(dy,2) + Math.pow(dz,2));
  return distance;
};

var coorDiff = function(coor1, coor2){
  var coor3 = [coor2[0]-coor1[0], coor2[1]-coor1[1], coor2[2]-coor1[2]];
  return coor3;
};

var coorAdd = function(coor1, coor2){
  var coor3 = [coor2[0]+coor1[0], coor2[1]+coor1[1], coor2[2]+coor1[2]];
  return coor3;
};

// Output Functions ------------------------------------------------------------
var recordSummary = function(IPsuccessful){
	var xdeg = xAngle / Math.PI * 180;
	var ydeg = yAngle / Math.PI * 180;
	var zdeg = zAngle / Math.PI * 180;
	console.log(randomPoint[0], randomPoint[1], randomPoint[2], xdeg, ydeg, zdeg, trialCount, structureCount, IPsuccessful);
};

var recordStructure = function(S1, S2, trialCount, structureCount){
	var significantFigure = 3;
	var decimal = Math.pow(10,significantFigure);
	console.log(S1.length + S2.length);
	console.log("IP T:" + trialCount + " S: " + structureCount);
	for(var i = 0; i < S1.length; i++){
		console.log("C "+ S1[i][0] + " " + S1[i][1] + " " + S1[i][2]);
	};
	for(var i = 0; i < S2.length; i++){
		x = Math.round(S2[i][0]*decimal)/decimal;
		y = Math.round(S2[i][1]*decimal)/decimal;
		z = Math.round(S2[i][2]*decimal)/decimal;
		console.log("O "+ x + " " + y + " " + z);
	};
	console.log("--------------------------------------------------------------");
};
// Quaternion Functions --------------------------------------------------------
function quaternion(w,x,y,z) {
  this.w = w;
  this.x = x;
  this.y = y;
  this.z = z;
  this.mult = function(Q1,Q2){
    var Q3 = new quaternion();
    Q3.w = Q1.w*Q2.w - Q1.x*Q2.x - Q1.y*Q2.y - Q1.z*Q2.z;
    Q3.x = Q1.x*Q2.w + Q1.w*Q2.x - Q1.z*Q2.y + Q1.y*Q2.z;
    Q3.y = Q1.y*Q2.w + Q1.z*Q2.x + Q1.w*Q2.y - Q1.x*Q2.z;
    Q3.z = Q1.z*Q2.w - Q1.y*Q2.x + Q1.x*Q2.y + Q1.w*Q2.z;
    return Q3;
  };
  this.conj = function(){
    var Qcon = new quaternion();
    Qcon.w = this.w;
    Qcon.x = -this.x;
    Qcon.y = -this.y;
    Qcon.z = -this.z;
    return Qcon;
  };
  this.inv = function(){
    var norm = Math.pow(this.w,2) + Math.pow(this.x,2) + Math.pow(this.y,2) + Math.pow(this.z,2);
    var Qinv = new quaternion();
    Qinv.w = this.w / norm;
    Qinv.x = -this.x / norm;
    Qinv.y = -this.y / norm;
    Qinv.z = -this.z / norm;
    return Qinv;
  };
  this.rotation = function(point,axisPoint1,axisPoint2,rotationAngle){
    var i = axisPoint2[0] - axisPoint1[0];
    var j = axisPoint2[1] - axisPoint1[1];
    var k = axisPoint2[2] - axisPoint1[2];
    var length = Math.sqrt( Math.pow(i,2) + Math.pow(j,2) + Math.pow(k,2) ); // What if length is 0??
    i = i / length;
    j = j / length;
    k = k / length; //console.log('ijk', i, j, k, "length:", length);
    var Qp = new quaternion();
    Qp.w = 0;
    Qp.x = point[0] - axisPoint2[0];
    Qp.y = point[1] - axisPoint2[1];
    Qp.z = point[2] - axisPoint2[2]; //console.log('Qp', Qp);
    var Qrot = new quaternion();
    Qrot.w = Math.cos(rotationAngle / 2.0);
    Qrot.x = Math.sin(rotationAngle / 2.0) * i;
    Qrot.y = Math.sin(rotationAngle / 2.0) * j;
    Qrot.z = Math.sin(rotationAngle / 2.0) * k; //console.log('Qrot', Qrot);
    var Q = this.mult(this.mult(Qrot,Qp),Qrot.inv()); //console.log('Q', Q);
    Q.x = Q.x + axisPoint2[0];
    Q.y = Q.y + axisPoint2[1];
    Q.z = Q.z + axisPoint2[2];
    return Q;
  };
};
// -----------------------------------------------------------------------------
// Define global variables
// -----------------------------------------------------------------------------
var initialCoordinates = [];  // Array of initial coordinates chosen from energy map if energy < (H)(H)
var energyMap = [];           // Array for storing energy values of each grid point in the unit cell
var newStructure = [];        // Array for storing coordinates of interpenetrating structure
var randomPoint;              // initial point for first atom selected randomly from initialCoordinates
var idx = 0;                  // Index for going through atoms in MOF unit cell
var structureCount = 0;       // Count the number of IP structures found
var trialCount = 0;           // Number of trials performed for interpenetration
var initialCoorTrialCount = 0;// To count the number of trials for different initial coordinates
var initialCoorIndex = 0;     // Index of coordinates selected from the energy map as inital coordinates
// -----------------------------------------------------------------------------
// ------------------------- Simulation Parameters -----------------------------
var cutOff = 13;              // Cut-off radius for LJ potential
var cubeLength = 26;          // Length of surrounding cube
var rotationFreedom = 12;			// Rotational freedom of interpenetrating structure (Ex: 60 degree rotation -> 360 / 6)
var rotationLimit = 20;

var energyLimitType = 'identical';
var energyScale = 1E10;				// Energy scale for collisionLimit and initialCoorEnergyLimit
var radiusScale = [0.25, 0.7, 0.6, 1.35];					// Emprical radii
var collisionLimit = []; var initialCoorEnergyLimit;
for(var i = 0; i < S1.length; i++){
	switch(energyLimitType){
		case 'sigmaScaled':
			collisionLimit[i] = (8.53/S1[i]) * energyScale;					  // Sigma Scaled
			initialCoorEnergyLimit = S1[1] * energyScale;
			break;
		case 'radiusScaled':
			collisionLimit[i] = (1.53/radiusScale[i]) * energyScale;		// Radius Scaled
			initialCoorEnergyLimit = radiusScale[3] * energyScale;
			break;
		case 'identical':
			collisionLimit[i] = 3*energyScale;									// Identical
			initialCoorEnergyLimit = 3 * energyScale;
			break;
	};
};
// ------------------------ Energy Map Variables -------------------------------
var grid = {};
grid.scale = 1; // The scale must be 1, 10, 100 ...
grid.size = grid.scale * cubeLength;
grid.length = cubeLength / grid.size;
grid.decimal = grid.size / cubeLength;
var eMapIndex;

var minCoor = -13; var maxCoor = 13;
// --------------------------- Rotation Variables ------------------------------
var q = new quaternion(0,1,1,1);
var xAngle, yAngle, zAngle, newX, newY, newZ;
var translationVector;
// -----------------------------------------------------------------------------
//------------------------- INITIALIZATION -------------------------------------
console.log('--------------------------------------------------------------------------------');
console.log('Energy Limit Type: ' + energyLimitType + ' Initial Coordinate Energy Limit: ' + initialCoorEnergyLimit);
console.log('Energy Limits: ' + 'H: ' + collisionLimit[0] + ' C: ' + collisionLimit[1] + ' O: ' + collisionLimit[2] + ' Zn: ' + collisionLimit[3]);
console.log('Rotational Freedom: ' + Math.round(360/rotationFreedom));
console.log('Rotation Limit: ' + rotationLimit);
console.log('Interpolation: ' + 'ON');
// Calculating the energy Map---------------------------------
calculateEnergyMap(difAtomsInfo, UC1, MOF5_Extended, cutOff);

console.log('Energy Map Calculated');
console.log('Total length: ' + energyMap.length + ' Grid size: ' + grid.length);
console.log('--------------------------------------------------------------------------------');

// Determine initial coordinates from energy map
for(var i = 0; i < energyMap.length; i++){
  if(energyMap[i][4] < initialCoorEnergyLimit){
    initialCoordinates.push([energyMap[i][0], energyMap[i][1], energyMap[i][2]]);
  };
};
// ------------------------ INTERPENETRATION -----------------------------------
var percentCompleted = 0;
var trialLimit = (initialCoordinates.length * rotationLimit);
var div = Math.round(trialLimit / 20);
var atomEnergy;
console.log('Interpenetration Starting with Trial Limit: ' + trialLimit);
console.log('--------------------------------------------------------------------------------');
for(var t = 0; t < trialLimit; t++){
	// Interpenetration trial loop
	loop1:
	for(var idx = 0; idx < MOF5.length; idx++){
	  if(idx === 0){
	    // For a number of trials for rotation it will start from the same initial point
	    if(trialCount % rotationLimit === 0){
	      randomPoint = initialCoordinates[initialCoorIndex];
	      initialCoorIndex++;
	      initialCoorTrialCount++;
	    };

	    // Select rotation axis and angle
	    xAngle = Math.PI*2*Math.floor(Math.random()*(rotationFreedom))/rotationFreedom;
	    yAngle = Math.PI*2*Math.floor(Math.random()*(rotationFreedom))/rotationFreedom;
	    zAngle = Math.PI*2*Math.floor(Math.random()*(rotationFreedom))/rotationFreedom;

			// Rotate first atom in the structure according to selected rotation
	    newCoor = [MOF5[0][0], MOF5[0][1], MOF5[0][2]];
	    q = q.rotation(newCoor, [0,0,0], [1,0,0], xAngle);
	    newCoor = [q.x, q.y, q.z];
	    q = q.rotation(newCoor, [0,0,0], [0,1,0], yAngle);
	    newCoor = [q.x, q.y, q.z];
	    q = q.rotation(newCoor, [0,0,0], [0,0,1], zAngle);
	    newCoor = [q.x, q.y, q.z];
	    trialCount++;

			// Calculate translation vector to translate each atom to the new rotated position
	    translationVector = coorDiff(newCoor, randomPoint);

			// Record coordinates for the new interpenetrating stucture
	    newStructure[structureCount] = [];
	    newStructure[structureCount][idx] = [randomPoint[0], randomPoint[1], randomPoint[2]];
	  };

		// Rotate the structure according to the selected rotation
    newCoor = [MOF5[idx][0], MOF5[idx][1], MOF5[idx][2]];
    q = q.rotation(newCoor, [0,0,0], [1,0,0], xAngle);
    newCoor = [q.x, q.y, q.z];
    q = q.rotation(newCoor, [0,0,0], [0,1,0], yAngle);
    newCoor = [q.x, q.y, q.z];
    q = q.rotation(newCoor, [0,0,0], [0,0,1], zAngle);
    newCoor = [q.x, q.y, q.z];

		// Translate the rotated structrue
    newCoor = coorAdd(newCoor, translationVector);

    // Apply periodic boundary conditions to find new coordinate in the unit cell (for energy map)
    pbcCoor = PBC(newCoor, 13);

		// Calculate index of the point in the energyMap array
    //eMapIndex = findEmapIndex(pbcCoor, grid.decimal, minCoor, maxCoor);

		// Find the type of atom being added
    atomType = findAtomType(MOF5[idx][3]);

		// 3D Interpolation for energy values
		atomEnergy = trInterpolate(pbcCoor, 4);

		// Check for energy level (collision)
    if(atomEnergy[atomType-3] >= collisionLimit[atomType-3]){
			//recordSummary(0);
			break loop1;
    } else{
      newStructure[structureCount][idx] = [newCoor[0], newCoor[1], newCoor[2]];
    };
		if(idx === MOF5.length-1){
	    //recordStructure(MOF5, newStructure[structureCount], trialCount, structureCount);
	    structureCount++;
			//recordSummary(1);
		};
	};
	if(t % div === 0){
		percentCompleted = Math.round(t / trialLimit * 100);
		console.log('Percent: ' + percentCompleted + ' Structure Count: ' + structureCount + ' Trial Count: ' + trialCount + ' IC: ' + initialCoorTrialCount);
	};
};
console.log('--------------------------------------------------------------------------------');
console.log('COMPLETED' + ' Structure Count: ' + structureCount + ' Trial Count: ' + trialCount + ' IC: ' + initialCoorTrialCount);
