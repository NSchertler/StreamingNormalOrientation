#include "utils.h"

void getSegmentColor(int segment, float& r, float& g, float& b)
{
	float H, S, V;
	H = 0.0f + 1.0f * ((3 * segment)  % 13) / 13.0f; //0.0 == 1.0
	S = 1.0;
	V = 0.5f + 0.5f * (segment % 2);

	//Calculate color from hue
	float _r = std::max(0.0f, std::min(1.0f, std::abs(H * 6 - 3) - 1));
	float _g = std::max(0.0f, std::min(1.0f, 2 - std::abs(H * 6 - 2)));
	float _b = std::max(0.0f, std::min(1.0f, 2 - std::abs(H * 6 - 4)));

	//Introduce saturation and brightness
	r = ((_r - 1.0f) * S + 1.0f) * V;
	g = ((_g - 1.0f) * S + 1.0f) * V;
	b = ((_b - 1.0f) * S + 1.0f) * V;
}
