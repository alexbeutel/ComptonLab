#include <iostream>
#include <fstream>
#include <vector>

void runCs137() {
	CTplotCalibration("1300sec--40deg-with.TKA", "1300sec--40deg-without.TKA",
		1000.0, 16000.0, 5700.0, 4000.0, 7000.0);
}
