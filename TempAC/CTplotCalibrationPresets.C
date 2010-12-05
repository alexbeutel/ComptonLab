#include <iostream>
#include <fstream>
#include <vector>

void runCs137() {
	CTplotCalibration("cal-co57.TKA", "cal-bg.TKA",
		1000.0, 16000.0, 2000.0, 1500.0, 3000.0);
}