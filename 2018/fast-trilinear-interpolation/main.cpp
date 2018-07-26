#define _USE_MATH_DEFINES
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>

using namespace std;

struct Color {
	float r, g, b, a;
};

Color operator*(const Color& c, float v)
{
	return Color{c.r*v, c.g*v, c.b*v};
}

Color operator*(float v, const Color& c)
{
	return Color{ c.r*v, c.g*v, c.b*v };
}

Color operator+(const Color& c1, const Color& c2)
{
	return Color{ c1.r+c2.r, c1.g+c2.g, c1.b+c2.b };
}

Color operator-(const Color& c1, const Color& c2)
{
	return Color{ c1.r - c2.r, c1.g - c2.g, c1.b - c2.b };
}

ostream& operator<<(ostream& os, Color c)
{
	os << '(' << c.r << ',' << c.g << ',' << c.b << ')';
	return os;
}

/* precompute trilinear coefficients
* v: f(x0,y0,z0), f(x1,y0,z0), f(x0,y1,z0), f(x1,y1,z0)
*    f(x0,y0,z1), f(x1,y0,z1), f(x0,y1,z1), f(x1,y1,z1)
*/
template <typename T>
array<T,8> PreComp(float x0, float y0, float z0,
	float x1, float y1, float z1,				
	const array<T,8>& v)
{
	const float epsilon = 10e-6f;
	if (x1 - x0 < epsilon) {
		x1 = 1.f;
		x0 = 0.f;
	}
	if (y1 - y0 < epsilon) {
		y1 = 1.f;
		y0 = 0.f;
	}
	if (z1 - z0 < epsilon) {
		z1 = 1.f;
		z0 = 0.f;
	}

	float deno = (x1 - x0)*(y1 - y0)*(z1 - z0);
	float nume = 1.f / deno;

	T a1 = x1 * v[0] - x0 * v[1];
	T a2 = x1 * v[2] - x0 * v[3];
	T a3 = x1 * v[4] - x0 * v[5];
	T a4 = x1 * v[6] - x0 * v[7];
	T b1 = v[1] - v[0];
	T b2 = v[3] - v[2];
	T b3 = v[5] - v[4];
	T b4 = v[7] - v[6];
	
	T c5 = y1 * a1 - y0 * a2;
	T c6 = y1 * a3 - y0 * a4;
	T d5 = a2 - a1;
	T d6 = a4 - a3;
	T e5 = y1 * b1 - y0 * b2;
	T e6 = y1 * b3 - y0 * b4;
	T f5 = b2 - b1;
	T f6 = b4 - b3;
	
	T h1 = z1 * c5 - z0 * c6;
	T h2 = z1 * d5 - z0 * d6;
	T h3 = z1 * e5 - z0 * e6;
	T h4 = z1 * f5 - z0 * f6;
	T h5 = c6 - c5;
	T h6 = d6 - d5;
	T h7 = e6 - e5;
	T h8 = f6 - f5;

	array<T, 8> h = {h1*nume, h2*nume, h3*nume, h4*nume, 
		h5*nume, h6*nume, h7*nume, h8*nume};
	return h;
}

vector<Color> FastTrilinear(const vector<Color>& voxel,
	int width, int height, int depth, float ratio)
{
	int dst_w = int(width * ratio);
	int dst_h = int(height * ratio);
	int dst_d = int(depth * ratio);

	auto f = [=](int x, int y, int z) {
		int idx = (min(z, depth - 1)*height + min(y, height - 1))*width
			+ min(x, width - 1);
		return voxel[idx];
	};
	vector<array<Color, 8>> coefs(voxel.size());

	int idx = 0;
	for (int z = 0; z < depth; z++) {
		float z0 = float(z)*dst_d/ depth;
		float z1 = float(z+1)*dst_d/ depth;
		for (int y = 0; y < height; y++) {
			float y0 = float(y)*dst_h/ height;
			float y1 = float(y+1)*dst_h/ height;
			for (int x = 0; x < width; x++) {
				float x0 = float(x)*dst_w/ width;
				float x1 = float(x+1)*dst_w/ width;
				array<Color, 8> vals = {
					f(x, y, z),		f(x + 1, y, z),
					f(x, y + 1, z),	f(x + 1, y + 1, z),
					f(x, y, z + 1),	f(x + 1, y, z + 1),
					f(x, y + 1,z + 1),f(x + 1,y + 1,z + 1)
				};
				coefs[idx] = PreComp<Color>(x0, y0, z0, x1, y1, z1, vals);
				idx++;
			}
		}
	}
	vector<Color> dst(dst_w*dst_h*dst_d);
	for (int i = 0; i < width*height*depth; i++) {
		int zc = i / (width*height);
		int yc = i / width % height;
		int xc = i % width;
		array<Color, 8> h = coefs[i];
		int zst = zc * dst_d / depth;
		int zed = (zc+1) * dst_d / depth;
		int yst = yc * dst_h / height;
		int yed = (yc+1) * dst_h / height;
		int xst = xc * dst_w / width;
		int xed = (xc+1) * dst_w / width;
		for (int z = zst; z < zed; z++) {
			float zf = float(z);
			for (int y = yst; y < yed; y++) {
				float yf = float(y);
				int idx_base = dst_w * (z*dst_h + y);
				for (int x = xst; x < xed; x++) {
					float xf = float(x);
					Color val = zf * (yf*(h[7] * xf + h[5]) + h[6] * xf + h[4])
						+ yf * (h[3] * xf + h[1]) + h[2] * xf + h[0];
					int idx = idx_base + x;
					dst[idx] = val;
				}
			}
		}
	}
	return dst;
}

vector<Color> GenVolume(int w, int h, int d)
{
	const float pi = float(M_PI);
	float R = 0.7f / 2;
	float r = 0.2f / 2;
	vector<Color> volume(w*h*d, Color{});
	for (float u = 0; u < 1; u += 0.01f) {
		float x = 0.4f*cos(u*4*pi)+0.5f;
		float y = 0.4f*sin(u *4 * pi)+0.5f;
		float z = u;
		Color c = {x,y,z, z};
		int cx = min(max(int(x*w), 0), w - 1);
		int cy = min(max(int(y*h), 0), h - 1);
		int cz = min(max(int(z*d), 0), d - 1);
		volume[cz*w*h + cy*w + cx] = c;
	}
	return volume;
}

void WriteVolume(string file, vector<Color> vol, int w, int h, int d)
{
	ofstream ofs(file, ios::binary);

	for (int i = 0; i < w*h*d; i++) {
		int zc = i / (w*w);
		int yc = i / w % h;
		int xc = i % w;
		Color c = vol[i];
		unsigned char ur = (unsigned char)min(max(c.r * 255, 0.f), 255.f);
		unsigned char ug = (unsigned char)min(max(c.g * 255, 0.f), 255.f);
		unsigned char ub = (unsigned char)min(max(c.b * 255, 0.f), 255.f);
		unsigned char ua = max(max(ur, ug), ub);
		unsigned char color[4] = { ub,ug,ur, ua};
		ofs << float(xc) / w << '\t'
			<< float(yc) / h << '\t'
			<< float(zc) / d << '\t'
			<< *((int*)color) << endl;
	}
}

void test()
{
	int w = 20, h = 20, d = 20;
	auto vol = GenVolume(w, h, d);
	auto vol_interp = FastTrilinear(vol, w, h, d, 1.f);
	float var = 0;
	int N = w * h*d;
	float maxd = 0;
	for (int i = 0; i < N; i++) {
		auto c1 = vol[i];
		auto c2 = vol_interp[i];
		auto d = c1 - c2;
		float sqsum = d.r*d.r + d.g*d.g + d.b*d.b + d.a*d.a;
		var += sqsum;
		maxd = max(maxd, sqsum);
	}
	var /= N;
	cout << "variance: " << var << endl;
	cout << "max diff: " << maxd << endl;
}

int main()
{
	test();

	int w = 20, h = 20, d = 20;
	float ratio = 5;
	int dst_w = int(w * ratio);
	int dst_h = int(h * ratio);
	int dst_d = int(d * ratio);

	auto vol = GenVolume(w, h, d);
	auto vol_interp = FastTrilinear(vol, w, h, d, ratio);
	WriteVolume("origin.txt", vol, w, h, d);
	WriteVolume("interp.txt", vol_interp, dst_w, dst_h, dst_d);

	return 0;
}
