#pragma once

#include <Eigen/Eigen>

// Mixin vertex construction

struct VBase
{
	void ReadFromString(char** s) { }
	void Transform(const Eigen::Matrix3f& transform) { }
	void Move(const Eigen::Vector3f& t) { }
	static void WritePLYProperties(std::ostream& s) { }
};

template <typename TBase>
struct VPosition : public TBase
{
	Eigen::Vector3f position;

	void ReadFromString(char** s)
	{
		TBase::ReadFromString(s);

		for (int i = 0; i < 3; ++i)
			position[i] = strtof(*s, s);
	}

	void Transform(const Eigen::Matrix3f& transform) 
	{
		TBase::Transform(transform);

		position = transform * position;
	}

	void Move(const Eigen::Vector3f& t) 
	{
		TBase::Move(t);

		position += t;
	}

	static void WritePLYProperties(std::ostream& s)
	{
		TBase::WritePLYProperties(s);

		s   << "property float x" << std::endl
			<< "property float y" << std::endl
			<< "property float z" << std::endl;
	}
};

template <typename TBase>
struct VNormal : public TBase
{
	Eigen::Vector3f normal;

	void ReadFromString(char** s)
	{		
		TBase::ReadFromString(s);

		for (int i = 0; i < 3; ++i)
			normal[i] = strtof(*s, s);
	}

	void Transform(const Eigen::Matrix3f& transform)
	{
		TBase::Transform(transform);

		normal = transform * normal;
	}

	void Move(const Eigen::Vector3f& t)
	{
		TBase::Move(t);
	}

	static void WritePLYProperties(std::ostream& s)
	{
		TBase::WritePLYProperties(s);

		s   << "property float nx" << std::endl
			<< "property float ny" << std::endl
			<< "property float nz" << std::endl;
	}
};

template <typename TBase>
struct VColor : public TBase
{
	Eigen::Matrix<unsigned char, 3, 1> color;

	void ReadFromString(char** s)
	{
		TBase::ReadFromString(s);

		for (int i = 0; i < 3; ++i)
			color[i] = (unsigned char)strtol(*s, s, 10);
	}

	void Transform(const Eigen::Matrix3f& transform)
	{
		TBase::Transform(transform);
	}

	void Move(const Eigen::Vector3f& t)
	{
		TBase::Move(t);
	}

	static void WritePLYProperties(std::ostream& s)
	{
		TBase::WritePLYProperties(s);

		s   << "property uchar red" << std::endl
			<< "property uchar green" << std::endl
			<< "property uchar blue" << std::endl;
	}
};

typedef VPosition<VBase> VertexPosition;
typedef VNormal<VertexPosition> VertexPositionNormal;
typedef VColor<VNormal<VertexPosition>> VertexPositionNormalColor;

extern std::ostream& operator<<(std::ostream& os, const VertexPosition& v);

extern std::ostream& operator<<(std::ostream& os, const VertexPositionNormal& v);

extern std::ostream& operator<<(std::ostream& os, const VertexPositionNormalColor& v);