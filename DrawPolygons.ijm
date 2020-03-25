/*
* By Michael Frank 
* March 2020
* Provided as-is from a request on the ImageJ Forum
* https://forum.image.sc/t/quantification-of-the-distribution-of-vesicles-in-different-sections-of-a-cell-fluorescence-intensity/34300
*/

// If you install this macro, it will be mapped to the F2 Key
macro "Draw Concentric Quadrants [F2]" { 
	
	// Clean previous ROIs
	roiManager("Reset");

	waitForUser("Draw an N-sided polygon around the object, the direction must be clockwise! (left to right)");

	// Get the line end points
	getSelectionCoordinates(x,y);


	// Calcualtion of the areas and masspoints
	//----------------------------------------------------------------------------------------------

	/* Calculate the high and low point of the polygon */
	lowX = 0; highX = 0; lowY = 0; highY = 0;

	for(i = 0; i < x.length; i++) {
		if(lowX == 0 || lowX > x[i]){
			lowX = x[i];
		}

		if(highX == 0 || highX < x[i]){
			highX = x[i];
		}

		if(lowY == 0 || lowY > y[i]){
			lowY = y[i];
		}

		if(highY == 0 || highY < y[i]){
			highY = y[i];
		}
	}

	/* Determine a pseudo "Center of mass" */

	pdX = (lowX + highX) / 2;
	pdY = (lowY + highY) / 2;


	/* Draw the polygon lines and calculate the areas and masspoints*/

	areas = newArray(x.length);
	massPointsX = newArray(x.length);
	massPointsY = newArray(x.length);

	for(a = 0; a < x.length; a++) {
		second_index = escapeIndex(x.length, a);
		//drawPolygonLine(a, second_index, "Polygon " + a);
		areas[a] = crossProduct(x[a], y[a], x[second_index], y[second_index], pdX, pdY);
		massPointsX[a] = (x[a] + x[second_index] + pdX) / 3;
		massPointsY[a] = (y[a] + y[second_index] + pdY) / 3;
	}


	// Draw the mass-point of the area
	mainX = mainMassPoint(addNumberList(areas), massPointsX, areas);
	mainY = mainMassPoint(addNumberList(areas), massPointsY, areas);

	makeLine(mainX, mainY, mainX, mainY);
	Roi.setName("Masspoint");
	roiManager("add");


	// Calcualtion and drawing of the quadrants
	//----------------------------------------------------------------------------------------------------------------------

	sX1 = 0; sX2 = 0; sY1 = 0; sY2 = 0;
	currentIndex = 0;

	// calculate the 90 degree interception from the masspoint to the outer-polygon-line
	//-----------------------------------------------------------------------------------
	for(i = 0; i < x.length; i++) {
		second_index = escapeIndex(x.length, i);
		if((y[i] < mainY || y[second_index] < mainY) && x[i] < mainX && x[second_index] > mainX){
			sX1 = x[i]; sX2 = x[second_index]; sY1 = y[i]; sY2 = y[second_index];
			currentIndex = i + 1;
		}
	}

	// calculate and draw the first intersection line
	iX = mainX;
	iY = calcStartPointOfIntersection(sX1, sY1, sX2, sY2, mainX);

	firstIPX = iX;
	firstIPY = iY;

	//drawIntersectionLine(mainX, mainY, iX, iY, "Intersection line 1");

	// calculate and draw the intersection lines
	//-----------------------------------------------------------------------------------
	sum_areas = addNumberList(areas); // sum the part-areas of the figures for the whole polygon area

	for(u = 0; u < 6; u++) {
		// determination boolean and area-holder
		determineLine = true;
		currentArea = 0;
		polygon_array_X = newArray(mainX, iX);
		polygon_array_Y = newArray(mainY, iY);

		// iterate through the polygon points until the area threshold is reached
		while(determineLine) {
			if(((crossProduct(iX, iY, mainX, mainY, x[currentIndex], y[currentIndex])) + currentArea) >= (sum_areas / 6)) {
				tX = iX;

				// case: the line goes up on the x-axis
				if(iX < x[currentIndex]) {
					while(tX < x[currentIndex]) {
						tX += 0.001;
						tY = calcStartPointOfIntersection(iX, iY, x[currentIndex], y[currentIndex], tX);

						if((crossProduct(mainX, mainY, iX, iY, tX, tY)) + currentArea >= (sum_areas / 6)) {
								iX = escapeFirstIntersectionPoint(tX, firstIPX, u); polygon_array_X = Array.concat(polygon_array_X, iX);
								iY = escapeFirstIntersectionPoint(tY, firstIPY, u); polygon_array_Y = Array.concat(polygon_array_Y, iY);
								drawPolygonPart(polygon_array_X, polygon_array_Y, "Polygon " + u + 1);
								determineLine = false;
								break;
						} 
					}

				// case: the line goes down on the x-axis
				} else {
					while(tX > x[currentIndex]) {
						tX -= 0.001;
						tY = calcStartPointOfIntersection(iX, iY, x[currentIndex], y[currentIndex], tX);

						if((crossProduct(mainX, mainY, iX, iY, tX, tY)) + currentArea >= (sum_areas / 6)) {
								iX = escapeFirstIntersectionPoint(tX, firstIPX, u); polygon_array_X = Array.concat(polygon_array_X, iX);
								iY = escapeFirstIntersectionPoint(tY, firstIPY, u); polygon_array_Y = Array.concat(polygon_array_Y, iY);
								drawPolygonPart(polygon_array_X, polygon_array_Y, "Polygon " + u + 1);
								determineLine = false;
								break;
						} 
					}
				}
			// case the area between the current intersection- and index point is smaller the 1/6 of the complete area
			} else {
				currentArea += (crossProduct(iX, iY, mainX, mainY, x[currentIndex], y[currentIndex]));
				iX = x[currentIndex]; polygon_array_X = Array.concat(polygon_array_X, iX);
				iY = y[currentIndex]; polygon_array_Y = Array.concat(polygon_array_Y, iY);
				currentIndex = escapeIndex(x.length, currentIndex);
			}
		}
	}
}

// escapes the OutOfBoundArrayException
function escapeIndex(length, currentIndex) {
	if(length - 1 == currentIndex) {
		return 0;
	} else {
		return (currentIndex + 1);
	}
}

function escapeFirstIntersectionPoint(cIP, fIP, index) {
	if(index == 5) {
		return fIP;
	} else {
		return cIP;
	}
}

// calculates the y-coordinate of the intersection point of a two-point equation
function calcStartPointOfIntersection(x1, y1, x2, y2, iX) {
	slope = (y2 - y1) / (x2 - x1); // Calculate the slope from the two points

	// Substitute the slope for 'm' in the slope intercept form of the equation

	b = y1 - (slope * x1); // Substitute either point into the equation: y = mx + b | while m is the calculated slope
	iY = slope * iX + b;	  // Substitute b in the equation

	return iY;
}

// calculates the mass-point of the polygon by the part-areas/-mass-points
function mainMassPoint(fullArea, polygon_areas, polygon_masspoints) {
	result = 0;
	for(i = 0; i < polygon_areas.length; i++) {
		result += polygon_masspoints[i] * polygon_areas[i];
	}
	result *= (1 / fullArea);
	return result;
}

// calculates the cross product of three points
function crossProduct(aX, aY, bX, bY, cX, cY) {
	abX = bX - aX; 
	abY = bY - aY; 

	acX = cX - aX; 
	acY = cY - aY; 

	result = abX * acY - acX * abY;

	// absolute value of the area
	if(result < 0) {
		result = (-1) * result;
	}
	return result
}

function addNumberList(list) {
	result = 0;
	for(i = 0; i < list.length; i++) {
		result += list[i];
	}
	return result;
}

function drawPolygonLine(index, second_index, name){
	makeLine(x[index], y[index], x[second_index], y[second_index]);
	Roi.setName(name)
	roiManager("add");
}

function drawIntersectionLine(x1, y1, x2, y2, name){
	makeLine(x1, y1, x2, y2);
	Roi.setStrokeColor("yellow");
	Roi.setName(name);
	roiManager("add");
}

/* 

									CODE OF SHAME 
	(the ImageJ function "makePolygon" takes N-Integer as the point-parameters, the problem is
	that a Polygon in this script can have N-Points, therefore we have to set the points manually!

	I have not much experience with this specific scripting-language, so if you know how refactor
	the block, feel free to create a pull-request. :) )
																									*/
function drawPolygonPart(paXs, paYs, name) {
	if(paXs.length == 1) {
		makePolygon(paXs[0], paYs[0]);
	} else if (paXs.length == 2){
		makePolygon(paXs[0], paYs[0], paXs[1], paYs[1]);
	} else if (paXs.length == 3){
		makePolygon(paXs[0], paYs[0], paXs[1], paYs[1], paXs[2], paYs[2]);
	} else if (paXs.length == 4){
		makePolygon(paXs[0], paYs[0], paXs[1], paYs[1], paXs[2], paYs[2], paXs[3], paYs[3]);		
	} else if (paXs.length == 5){
		makePolygon(paXs[0], paYs[0], paXs[1], paYs[1], paXs[2], paYs[2], paXs[3], paYs[3], paXs[4], paYs[4]);		
	} else if (paXs.length == 6){
		makePolygon(paXs[0], paYs[0], paXs[1], paYs[1], paXs[2], paYs[2], paXs[3], paYs[3], paXs[4], paYs[4], paXs[5], paYs[5]);		
	} else if (paXs.length == 7){
		makePolygon(paXs[0], paYs[0], paXs[1], paYs[1], paXs[2], paYs[2], paXs[3], paYs[3], paXs[4], paYs[4], paXs[5], paYs[5], paXs[6], paYs[6]);		
	} else if (paXs.length == 8){
		makePolygon(paXs[0], paYs[0], paXs[1], paYs[1], paXs[2], paYs[2], paXs[3], paYs[3], paXs[4], paYs[4], paXs[5], paXs[5], paXs[6], paYs[6], paYs[7], paYs[7]);		
	} 
	Roi.setStrokeColor("yellow");
	Roi.setName(name);
	roiManager("add");
}
