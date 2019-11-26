#pragma once

#include <cmath>
#include <cassert>

/* Классы и функции пространственной геометрии */
namespace Geom3D {

/* Величина погрешности */
const double epsilon{ 0.00001 };

/* Квадрат числа */
template <typename T = double>
T sqr(T x) { return x * x; }

/* Точка */
template <typename T = double>
struct Point;

/* Направляющий вектор */
template <typename T = double>
struct Vector;

/* Точка */
template <typename T>
struct Point
{
	Point(const T x = 0, const T y = 0, const T z = 0)
		: x(x)
		, y(y)
		, z(z) 
	{}

	Point(const Point<T>& p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
	}

	/* Вычитание */
	Vector<T> operator - (const Point<T>& p) const
	{
		return Vector<T>(x - p.x, y - p.y, z - p.z);
	}

	/* Сложение с вектором */
	Point<T> operator + (const Vector<T>& v) const
	{
		return Point<T>(x + v.l, y + v.m, z + v.n);
	}

	/* Вычитание вектора */
	Point<T> operator - (const Vector<T>& v) const
	{
		return Point<T>(x - v.l, y - v.m, z - v.n);
	}

	/* Компоненты точки */
	T x;
	T y;
	T z;
};

/* Направляющий вектор */
template <typename T>
struct Vector
{
	Vector(const Point<T>& p0, const Point<T>& p1)
		: l(p1.x - p0.x)
		, m(p1.y - p0.y)
		, n(p1.z - p0.z)
	{}

	Vector(T x, T y, T z)
		: l(x)
		, m(y)
		, n(z)
	{}

	Vector(const Vector<T>& v)
	{
		l = v.l;
		m = v.m;
		n = v.n;
	}

	/* Сложение */
	Vector<T> operator + (const Vector<T>& v) const
	{
		return Vector<T>(l + v.l, m + v.m, n + v.n);
	}

	/* Вычитание */
	Vector<T> operator - (const Vector<T>& v) const
	{
		return Vector<T>(l - v.l, m - v.m, n - v.n);
	}

	/* Умножение на скаляр */
	Vector<T> operator * (T s) const
	{
		return Vector<T>(l * s, m * s, n * s);
	}

	/* Деление на скаляр */
	Vector<T> operator / (T s) const
	{
		return Vector<T>(l / s, m / s, n / s);
	}

	/* Умножение на вектор (векторное произведение) */
	Vector<T> operator * (const Vector<T>& v) const
	{
		return Vector<T>(m * v.n - n * v.m, n * v.l - l * v.n, l * v.m - m * v.l);
	}

	/* Скалярное произведение */
	T dot(const Vector<T>& v) const
	{
		return l * v.l + m * v.m + n * v.n;
	}

	/* Квадрат длины вектора */
	T sqrlen() const { return sqr(l) + sqr(m) + sqr(n); }

	/* Длина вектора */
	T norm() const { return sqrt(sqrlen()); }	

	/* Нормированный вектор */
	Vector<T> unit() const
	{
		return *this / norm();
	}

	/* Компоненты вектора */
	T l;
	T m;
	T n;
};

/* Скалярное произведение двух векторов */
template <typename T>
T dot(const Vector<T>& u, const Vector<T>& v)
{
	return u.l * v.l + u.m * v.m + u.n * v.n;
}

/* Норма вектора (длина) */
template <typename T>
T norm(const Vector<T>& v)
{
	return sqrt(dot(v, v));
}

/* Отрезок */
template <typename T = double>
class Segment
{
public:
	Segment(const Point<T>& p0, const Point<T>& p1)
		: m_p0(p0)
		, m_p1(p1)
		, m_v(p0, p1)
	{}

	Segment(const Point<T>& p, const Vector<T>& v)
		: m_p0(p)
		, m_v(v)
		, m_p1(p + v)
	{}

	Segment(const Segment<T>& s)
	{
		m_p0 = s.m_p0;
		m_p1 = s.m_p1;
		m_v = s.m_v;
	}

	/* Точка прямой, заданной параметрически:
	 * p(t) = p0 + v * t, где 0<t<1
	 */
	Point<T> point(T t) const
	{
		assert(t <= 1);
		assert(t >= 0);

		return m_p0 + m_v * t;
	}

	/* Направляющий вектор */
	const Vector<T>& vector() const
	{
		return m_v;
	}

	/* Признак параллельности с другим отрезком */
	bool parallel(const Segment<T>& s) const
	{
		return (m_v * s.m_v).norm() < epsilon;
	}

	/* Точка начала отрезка */
	const Point<T>& begPoint() const
	{
		return m_p0;
	}

	/* Точка конца отрезка */
	const Point<T>& endPoint() const
	{
		return m_p1;
	}

private: /* Кэшированные значения */

	/* Начальная точка */
	Point<T> m_p0;

	/* Концевая точка */
	Point<T> m_p1;

	/* Направляющий вектор */
	Vector<T> m_v;
};

} // namespace Geom