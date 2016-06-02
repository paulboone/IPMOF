// -------------------------------------------------------------------------------------------------
// Energy Map Related Functions --------------------------------------------------------------------
GLOBAL.findEmapIndex = function(coor, decimal, eMapMax, eMapMin){
	var sideLength = []; var rndCoor = [];
	for(var i = 0; i < coor.length; i++){
		sideLength.push(eMapMax[i] - eMapMin[i] + 1);
		rndCoor.push(Math.round(coor[i]*decimal)/decimal);
	};
	var eMapIndex = rndCoor[0]*sideLength[1]*sideLength[2] + rndCoor[1]*sideLength[2] + rndCoor[2];
	return eMapIndex;
};

GLOBAL.findAtomType = function(atomSymbol, eMapAtomNames, eMapAtomIndex){
	var atomIndex;
	for(var atom_idx = 0; atom_idx < eMapAtomNames.length; atom_idx++){
		if(atomSymbol === eMapAtomNames[atom_idx]){
			atomIndex = eMapAtomIndex[atom_idx];
		};
	};
	return atomIndex;
};

// Uses GLOBAL variabes of eMapAtomNames and eMapAtomIndex
GLOBAL.selectInitialCoordinates = function(referenceAtom, eMap, eLimit, fracUCV){
  var refAtomIndex = findAtomType(referenceAtom, eMapAtomNames, eMapAtomIndex);
  if(refAtomIndex === undefined){
    refAtomIndex = 3; console.log(referenceAtom, ' atom not found!');
  };
	var initialCoor = []; var eMapCoor, fracCoor, pbcCoor, pbcX, pbcY, pbcZ;
  for(var i = 0; i < eMap.length; i++){
		eMapCoor = [eMap[i][0], eMap[i][1], eMap[i][2]];
		fracCoor = car2frac(eMapCoor, baseMOF_UCsize, baseMOF_UCangle, fracUCV);
		pbcCoor = fracPBC(fracCoor);
		pbcCoor = frac2car(pbcCoor, baseMOF_UCsize, baseMOF_UCangle, fracUCV);
		pbcX = Math.round(eMapCoor[0]*10)/10;
		pbcY = Math.round(eMapCoor[1]*10)/10;
		pbcZ = Math.round(eMapCoor[2]*10)/10;
		if(pbcX === eMapCoor[0] && pbcY === eMapCoor[1] && pbcZ === eMapCoor[2]){
			if(eMap[i][refAtomIndex] < eLimit){
				initialCoor.push([eMap[i][0], eMap[i][1], eMap[i][2]]);
			};
		};
	};
	return initialCoor;
};

// PBC and Interpolation Functions -----------------------------------------------------------------
// 3D Linear Interpolation
// MODIFIES INPUT COORDINATE IF A ROUNDED COORDINATE VALUE IS GIVEN!!!!!!!!!!!!!
GLOBAL.trInterpolate = function(coor, atomIndex, eMap, eMapMax, eMapMin, gridSize){
	var coor1 = []; var coor0 = []; var dif = [];

	for(var i = 0; i < coor.length; i++){
		if(Math.round(coor[i]*10)/10 === coor[i]){
			coor[i] += 1E-10;
		};
		coor0[i] = Math.floor(coor[i]);
		coor1[i] = Math.ceil(coor[i]);
		dif[i] = ( coor[i] - coor0[i] ) / ( coor1[i] - coor0[i] );
	};

	var i000 = findEmapIndex(coor0, gridSize, eMapMax, eMapMin);
	var i100 = findEmapIndex([coor1[0],coor0[1],coor0[2]], gridSize, eMapMax, eMapMin);
	var i001 = findEmapIndex([coor0[0],coor0[1],coor1[2]], gridSize, eMapMax, eMapMin);
	var i101 = findEmapIndex([coor1[0],coor0[1],coor1[2]], gridSize, eMapMax, eMapMin);
	var i010 = findEmapIndex([coor0[0],coor1[1],coor0[2]], gridSize, eMapMax, eMapMin);
	var i110 = findEmapIndex([coor1[0],coor1[1],coor0[2]], gridSize, eMapMax, eMapMin);
	var i011 = findEmapIndex([coor0[0],coor1[1],coor1[2]], gridSize, eMapMax, eMapMin);
	var i111 = findEmapIndex(coor1, gridSize, eMapMax, eMapMin);

	var c00 = eMap[i000][atomIndex]*(1-dif[0]) + eMap[i100][atomIndex]*dif[0];
	var c01 = eMap[i001][atomIndex]*(1-dif[0]) + eMap[i101][atomIndex]*dif[0];
	var c10 = eMap[i010][atomIndex]*(1-dif[0]) + eMap[i110][atomIndex]*dif[0];
	var c11 = eMap[i011][atomIndex]*(1-dif[0]) + eMap[i111][atomIndex]*dif[0];

	var c0 = c00*(1-dif[1]) + c10*dif[1];
	var c1 = c01*(1-dif[1]) + c11*dif[1];

	var c = c0*(1-dif[2]) + c1*dif[2];

	return c;
};

GLOBAL.fracVolume = function(UCangle){
	var alp = degToRad(UCangle[0]);
	var bet = degToRad(UCangle[1]);
	var gam = degToRad(UCangle[2]);

	var v = 1 - Math.pow(Math.cos(alp),2) - Math.pow(Math.cos(bet),2);
	v += - Math.pow(Math.cos(gam),2) + 2*Math.cos(alp)*Math.cos(bet)*Math.cos(gam);
	v = Math.sqrt(v);

	return v;
};


// Output Functions --------------------------------------------------------------------------------
GLOBAL.recordSummary = function(IPsuccessful){
	console.log(firstPoint[0], firstPoint[1], firstPoint[2], xAngle, yAngle, zAngle, trialCount, structureCount, IPsuccessful);
};

GLOBAL.recordStructure = function(S1, S2, trialCount, structureCount, xAngle, yAngle, zAngle, structureTotalEnergy){
	var significantFigure = 3;
	var decimal = Math.pow(10,significantFigure);
	var xdeg = Math.round(xAngle / Math.PI * 180);
	var ydeg = Math.round(yAngle / Math.PI * 180);
	var zdeg = Math.round(zAngle / Math.PI * 180);
	console.log(S1.length + S2.length);
	console.log("Structure: " + structureCount + " IP Trial: " + trialCount);
	console.log("Rotation x: " + xdeg + " y: " + ydeg + " z: " + zdeg);
	console.log("Energy: " + structureTotalEnergy);
	var x, y, z;
	for(var i = 0; i < S1.length; i++){
		x = Math.round(S1[i][0]*decimal)/decimal;
		y = Math.round(S1[i][1]*decimal)/decimal;
		z = Math.round(S1[i][2]*decimal)/decimal;
		console.log(S1[i][3] + " " + x + " " + y + " " + z);
	};
	var x2, y2, z2;
	for(var i = 0; i < S2.length; i++){
		x2 = Math.round(S2[i][0]*decimal)/decimal;
		y2 = Math.round(S2[i][1]*decimal)/decimal;
		z2 = Math.round(S2[i][2]*decimal)/decimal;
		console.log(S2[i][3] + " " + x2 + " " + y2 + " " + z2);
	};
	console.log("--------------------------------------------------------------");
};
