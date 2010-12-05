#include <iostream>
#include <fstream>
#include <vector>

void runCs137() {
	CTplotCalibration("60sec-0deg-with.TKA", "60sec-0deg-without.TKA",
		1000.0, 16000.0, 7000.0, 6000.0, 9000.0);
}