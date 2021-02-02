#pragma once
#include <iostream>


void HardThreshold(float &value,float threshold)
{
	if (abs(value)<threshold)
	{
		value = 0;
	}
}