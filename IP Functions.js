// ----------------------- Animation Functions ---------------------------------
var black = 0x000000;
var blue = 0x0055ff;
var red = 0xcc0000;
var green = 0x009933;
var purple = 0x990099;
var lineWidth = 100;

var addAtom = function(coor, atomRadius, atomColor){
	var atomGeo = new THREE.SphereGeometry(atomRadius, 8, 8);
	var atomMat = new THREE.MeshLambertMaterial( { color: atomColor } );
	newAtom = new THREE.Mesh(atomGeo, atomMat);
	newAtom.position.set( coor[0], coor[1], coor[2] );
	scene.add( newAtom );
};

var addLine = function(coor1, coor2, lineColor){
        var lineMat = new THREE.LineBasicMaterial({color: lineColor, linewidth: lineWidth});
        var lineGeo = new THREE.Geometry();

        lineGeo.vertices.push(
	        new THREE.Vector3(coor1[0], coor1[1], coor1[2]),
	        new THREE.Vector3(coor2[0], coor2[1], coor2[2])
        );

        var line = new THREE.Line( lineGeo, lineMat );
        scene.add( line );
};
// -----------------------------------------------------------------------------
var reset = function(objectNumber){
	idx = 0;
	for( var i = scene.children.length - 1; i > objectNumber; i--) {
		obj = scene.children[i];
		scene.remove(obj);
	};
};

var degToRad = function(deg){
	return deg/180*Math.PI;
};

var radToDeg = function(rad){
	return rad*180/Math.PI;
};

var findEmapIndex = function(coor, decimal, minCoor, maxCoor){
	var gridSize = Math.abs(minCoor) + Math.abs(maxCoor);
	for(var i = 0; i < 3; i++){
		coor[i] = Math.round(coor[i]*decimal)/decimal;
	};
	var eMapIndex = (coor[0]-minCoor)*Math.pow(gridSize+1,2) + (coor[1]-minCoor)*(gridSize+1) + coor[2]-minCoor;
	return eMapIndex;
};

var PBC = function(coor, cutOff){
	var pbcCoor = [];
	for(var i = 0; i < coor.length; i++){
		if(Math.abs(coor[i]) > cutOff){
			pbcCoor[i] = coor[i] - Math.sign(coor[i])*(2*cutOff);
		} else{
			pbcCoor[i] = coor[i];
		}
	};
	return pbcCoor;
};

var PBC2 = function(coor, cutOff){
	var pbcCoor = [];
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
	console.log(randomPoint[0], randomPoint[1], randomPoint[2], xAngle, yAngle, zAngle, trialCount, structureCount, IPsuccessful);
};

var recordStructure = function(S1, S2, trialCount, structureCount){
	significantFigure = 3;
	decimal = Math.pow(10,significantFigure);
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
var calculateEnergyMap = function(extendedStructure, ucA, ucB, ucC, cutOff){
	var E = [22.156, 52.873, 30.213, 62.441]; // Epsilon for H, C, O, Zn
	var S = [2.571, 3.431, 3.118, 2.462];     // Sigma for H, C, O, Zn
	var A = ['H', 'C', 'O', 'Zn'];
	var V = [];
	var Vtotal = [];
	var i = 0;
	minA = Math.round(- ucA / 2); maxA = Math.round(ucA / 2);
	minB = Math.round(- ucB / 2); maxB = Math.round(ucB / 2);
	minC = Math.round(- ucC / 2); maxC = Math.round(ucC / 2);

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
//------------------------------------------------------------------------------

  // Trilinear Interpolation
var trInterpolate = function(coor, numAtoms){
  var coor1 = []; var coor0 = []; var dif = [];

  for(var i = 0; i < coor.length; i++){
    coor0[i] = Math.floor(coor[i]);
    coor1[i] = Math.ceil(coor[i]);
    dif[i] = ( coor[i] - coor0[i] ) / ( coor1[i] - coor0[i] );
  };

  i000 = findEmapIndex(coor0, grid.decimal, minCoor, maxCoor);
  i100 = findEmapIndex([coor1[0],coor0[1],coor0[2]], grid.decimal, minCoor, maxCoor);
  i001 = findEmapIndex([coor0[0],coor0[1],coor1[2]], grid.decimal, minCoor, maxCoor);
  i101 = findEmapIndex([coor1[0],coor0[1],coor1[2]], grid.decimal, minCoor, maxCoor);
  i010 = findEmapIndex([coor0[0],coor1[1],coor0[2]], grid.decimal, minCoor, maxCoor);
  i110 = findEmapIndex([coor1[0],coor1[1],coor0[2]], grid.decimal, minCoor, maxCoor);
  i011 = findEmapIndex([coor0[0],coor1[1],coor1[2]], grid.decimal, minCoor, maxCoor);
  i111 = findEmapIndex(coor1, grid.decimal, minCoor, maxCoor);

  var V = [];
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
}; // 3D Interpolation Function
