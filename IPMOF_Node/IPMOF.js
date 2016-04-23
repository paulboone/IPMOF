// ------------------------------------ Simulation Parameters --------------------------------------
// Grid size of the energy map (Default value is 1)
var gridSize = 1;
// Rotation degree is select according to this parameter
// Possible rotation degrees are calculated as 360 / rotation freedom
var rotationFreedom = 6;
// Energy limit that is used while checking for collision during interpenetration
// Energy scale is used to scale the avg energy limit calculated from eMap (Default: 0.1)
var energyScale = 0.1;
var energyLimit = avgEnergyLimit * energyScale;
// Energy limit used in eliminating 'bad' points in the energy map
var initialCoorEnergyLimit = avgEnergyLimit;
// Number of random rotations performed for each trial point
var rotationLimit = 30;
// Outputing structure coordinates (true or false)
var outputStructures = true;
// ---------------------------Interpenetration Global Variables ------------------------------------
var structureTotalEnergy = 0; // Total potential energy of the interpenetrating MOF
var newStructure = [];        // Array for storing coordinates of interpenetrating structure
var firstPoint;               // Initial point for first atom selected randomly from initialCoordinates
var idx = 0;                  // Index for going through atoms in MOF unit cell
var structureCount = 0;       // Count the number of IP structures found
var trialCount = 0;           // Number of trials performed for interpenetration
var initialCoorTrialCount = 0;// To count the number of trials for different initial coordinates
var initialCoorIndex = 0;     // Index of coordinates selected from the energy map as inital coordinateS
var eMapIndex;                // Index of atom in the energy map
// Lower and upper coordinate limits of the energy map - used for finding energy map index
var eMapMax = [eMap[eMap.length-1][0], eMap[eMap.length-1][1], eMap[eMap.length-1][2]];
var eMapMin = [eMap[0][0], eMap[0][1], eMap[0][2]];
// Calculate fractional unit cell volume
var fracUCV = fracVolume(baseMOF_UCangle);
// Determine initial coordinates from energy map by omitting points with e > initialCoorEnergyLimit
// and omitting points that are outside the unit cell (checked by applying PBC)
var initialCoordinates = selectInitialCoordinates('C', eMap, initialCoorEnergyLimit, fracUCV);
// Quaternion Functions ----------------------------------------------------------------------------
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
    k = k / length;
    var Qp = new quaternion();
    Qp.w = 0;
    Qp.x = point[0] - axisPoint2[0];
    Qp.y = point[1] - axisPoint2[1];
    Qp.z = point[2] - axisPoint2[2];
    var Qrot = new quaternion();
    Qrot.w = Math.cos(rotationAngle / 2.0);
    Qrot.x = Math.sin(rotationAngle / 2.0) * i;
    Qrot.y = Math.sin(rotationAngle / 2.0) * j;
    Qrot.z = Math.sin(rotationAngle / 2.0) * k;
    var Q = this.mult(this.mult(Qrot,Qp),Qrot.inv());
    Q.x = Q.x + axisPoint2[0];
    Q.y = Q.y + axisPoint2[1];
    Q.z = Q.z + axisPoint2[2];
    return Q;
  };
};
// ------------------------------------ Rotation Variables -----------------------------------------
var q = new quaternion(0,1,1,1);
var xAngle, yAngle, zAngle, newX, newY, newZ;
var translationVector;
// -------------------------------------------------------------------------------------------------
// ************************************ INITIALIZATION *********************************************
var simulationSummary = "";
simulationSummary += '--------------------------------------------------------------------------------' + '\n';
simulationSummary += 'Base MOF Atoms: ' + baseMOF.length + ' Atom Names: ' + baseMOFatomNames + '\n';
simulationSummary += 'Base MOF Unit Cell a: ' + baseMOF_UCsize[0] + ' b: ' + baseMOF_UCsize[1] + ' c: ' + baseMOF_UCsize[2] + '\n';
simulationSummary += 'Base MOF Unit Cell alpha: ' + baseMOF_UCangle[0] + ' beta: ' + baseMOF_UCangle[1] + ' gamma: ' + baseMOF_UCangle[2] + '\n';
simulationSummary += 'Mobile MOF Atoms: ' + MOF.length + ' Atom Names: ' + MOFatomNames + '\n';
simulationSummary += 'Mobile MOF Unit Cell a: ' + MOF_UCsize[0] + ' b: ' + MOF_UCsize[1] + ' c: ' + MOF_UCsize[2] + '\n';
simulationSummary += 'Mobile MOF Unit Cell alpha: ' + MOF_UCangle[0] + ' beta: ' + MOF_UCangle[1] + ' gamma: ' + MOF_UCangle[2] + '\n';
simulationSummary += 'Avg Energy Limit (eMap): ' + avgEnergyLimit + '\n';
simulationSummary += 'Energy Scale: ' + energyScale + '\n';
simulationSummary += 'Initial Coordinate Energy Limit: ' + initialCoorEnergyLimit + '\n';
simulationSummary += 'Energy Limits: ' + energyLimit + '\n';
simulationSummary += 'Rotational Freedom: ' + Math.round(360/rotationFreedom) + '\n';
simulationSummary += 'Rotation Limit: ' + rotationLimit + '\n';
simulationSummary += '--------------------------------------------------------------------------------' + '\n';
// *********************************** INTERPENETRATION ********************************************
var percentCompleted = 0;
var omittedCoordinates = eMap.length - initialCoordinates.length;
var trialLimit = (initialCoordinates.length * rotationLimit);
var div = Math.round(trialLimit / 19);
simulationSummary += 'Number of omitted candidate coordinates: ' + omittedCoordinates + '\n';
simulationSummary += 'Interpenetration Starting with Trial Limit: ' + trialLimit + '\n';
simulationSummary += '--------------------------------------------------------------------------------' + '\n';
for(var t = 0; t < trialLimit; t++){
  // Interpenetration trial loop
	loop1:
	for(var idx = 0; idx < MOF.length; idx++){
    if(idx === 0){

      if(trialCount % rotationLimit === 0){
        firstPoint = initialCoordinates[initialCoorIndex];
        initialCoorIndex++;
        initialCoorTrialCount++;
      };

      // Select rotation axis and angle
      xAngle = Math.PI*2*Math.floor(Math.random()*(rotationFreedom))/rotationFreedom;
      yAngle = Math.PI*2*Math.floor(Math.random()*(rotationFreedom))/rotationFreedom;
      zAngle = Math.PI*2*Math.floor(Math.random()*(rotationFreedom))/rotationFreedom;

      // Rotate first atom of the mobile MOF
      atomName = MOF[idx][3];
      newCoor = [MOF[idx][0], MOF[idx][1], MOF[idx][2]];
      q = q.rotation(newCoor, [0,0,0], [1,0,0], xAngle);
      newCoor = [q.x, q.y, q.z];
      q = q.rotation(newCoor, [0,0,0], [0,1,0], yAngle);
      newCoor = [q.x, q.y, q.z];
      q = q.rotation(newCoor, [0,0,0], [0,0,1], zAngle);
      newCoor = [q.x, q.y, q.z];
      trialCount++;

      translationVector = coorDiff(newCoor, firstPoint);

      newStructure[structureCount] = [];
      newStructure[structureCount][idx] = [firstPoint[0], firstPoint[1], firstPoint[2], atomName];
      idx++;
    };

    if(idx < MOF.length && idx > 0){
      atomName = MOF[idx][3];
      newCoor = [MOF[idx][0], MOF[idx][1], MOF[idx][2]];
      q = q.rotation(newCoor, [0,0,0], [1,0,0], xAngle);
      newCoor = [q.x, q.y, q.z];
      q = q.rotation(newCoor, [0,0,0], [0,1,0], yAngle);
      newCoor = [q.x, q.y, q.z];
      q = q.rotation(newCoor, [0,0,0], [0,0,1], zAngle);
      newCoor = [q.x, q.y, q.z];

      newCoor = coorAdd(newCoor, translationVector);
      fracCoor = car2frac(newCoor, baseMOF_UCsize, baseMOF_UCangle, fracUCV);
      pbcCoor = fracPBC(fracCoor);

      pbcCoor = frac2car(pbcCoor, baseMOF_UCsize, baseMOF_UCangle, fracUCV);
      eMapIndex = findEmapIndex(pbcCoor, gridSize, eMapMax, eMapMin);

      atomIndex = findAtomType(MOF[idx][3], eMapAtomNames, eMapAtomIndex);
      pointEnergy = trInterpolate(pbcCoor, atomIndex, eMap, eMapMax, eMapMin, gridSize);
      structureTotalEnergy += pointEnergy;
      if(pointEnergy > energyLimit){
        structureTotalEnergy = 0;
        break loop1;
      } else{
        newStructure[structureCount][idx] = [newCoor[0], newCoor[1], newCoor[2], atomName];
      };
    };

    if(idx === MOF.length-1) {
      structureTotalEnergy = Math.round(structureTotalEnergy/1E4)*1E4;
      structureTotalEnergy = structureTotalEnergy.toExponential();
      if(outputStructures){
        recordStructure(baseMOF, newStructure[structureCount], trialCount, structureCount+1, xAngle, yAngle, zAngle, structureTotalEnergy);
      };
      structureCount++;
      structureTotalEnergy = 0;
    };
  };

  if(t % div === 0){
    // Record simulation summary line
    percentCompleted = Math.round(t / trialLimit * 100);
    var tcText = trialCount;
    var totalDigit = trialLimit.toString().length - trialCount.toString().length;
    for(var d = 0; d < totalDigit; d++){ tcText += ' '; };
    rotationTrial = Math.round(t/rotationLimit) % rotationLimit;
    simulationSummary += 'Percent:  ' + percentCompleted + ' \tStructure Count: ' + structureCount;
    simulationSummary += ' \tTrial Count: ' + tcText + ' \tIC Count: ' + initialCoorTrialCount + '\n';
  };
};
// *********************************** Export Simulation Summary ***********************************
simulationSummary += 'Percent: 100' + '\tStructure Count: ' + structureCount;
simulationSummary += '\tTrial Count: ' + trialCount + ' \tIC Count: ' + initialCoorTrialCount + '\n';
simulationSummary += '--------------------------------------------------------------------------------' + '\n';
console.log(simulationSummary);
