#include <algorithm>
#include <limits>

#include "geom.h"
#include "print.h"

using namespace Geom3D;

/* Нахождение расстояния между двумя отрезками */
template <typename T>
T distance(const Segment<T>& s1, const Segment<T>& s2)
{
	/* Кэширование данных */
	auto uu = dot(s1.vector(), s1.vector());
	auto uv = dot(s1.vector(), s2.vector());
	auto vv = dot(s2.vector(), s2.vector());
	auto uw = dot(s1.vector(), s1.begPoint() - s2.begPoint());
	auto vw = dot(s2.vector(), s1.begPoint() - s2.begPoint());

	auto denominator = (uv * uv - uu * vv);

	/* Параметры, определяющие расположения ближайших точек на отрезках,
	 * заданных параметрически:
	 * Pc = P0 + sc * u
	 * Qc = Q0 + tc * v
	 */
	T sc = 0.0;
	T tc = 0.0;

	/* Решение */
	if (s1.parallel(s2))
	{
		/* Отрезки параллельны (особый случай) */
		if (uw < 0)
		{
			tc = 0.0;
			sc = -uw / uu;
		}
		else 
		{
			sc = 0.0;
			tc = uw / uv;
		}
	}
	else
	{
		/* Отрезки не параллельны */
		sc = (vv * uw - uv * vw) / denominator;
		tc = (uv * uw - uu * vw) / denominator;
	}

	/* Если ближайшие точки лежат внутри отрезков (отрезки перекрываются) */
	if (sc >= 0.0 && sc <= 1.0 && tc >= 0.0 && tc <= 1.0)
	{
		/* Искомое расстояние */
		return norm(s1.point(sc) - s2.point(tc));
	}
	else
	{
		/* Ищем ближайшие точки на концах отрезков */
		bool perpendicular = false;

		/* Если точка первого отрезка вышла за границы, подводим её к ближайшему краю */
		if (sc < 0)
		{
			sc = 0.0;

			/* Пересчитываем вторую точку, если она попадает в диапазон */
			if (uw > 0 && uw < uv)
			{
				tc = uw / uv;
				perpendicular = true;
			}
		}
		else if (sc > 1)
		{
			sc = 1.0;

			/* Пересчитываем вторую точку, если она попадает в диапазон */
			if ((uu > -uw) && (uu + uw < uv))
			{
				tc = (uu + uw) / uv;
				perpendicular = true;
			}
		}

		/* Если точка второго отрезка вышла за границы, подводим её к ближайшему краю */
		if (tc < 0)
		{
			tc = 0;

			/* Пересчитываем первую точку, если она попадает в диапазон */
			if ((-uw > 0) && (-uw < uu))
			{
				sc = -uw / uu;
				perpendicular = true;
			}
		}
		else if (tc > 1)
		{
			tc = 1;

			/* Пересчитываем первую точку, если она попадает в диапазон */
			if ((vv > vw) && (vv - vw < uv));
			{
				sc = (vv - vw) / uv;
				perpendicular = true;
			}
		}

		/* Если ни из одной из крайних точек невозможно построить перпендикуляр */
		if (!perpendicular)
		{
			/* Искомое расстояние: минимальное расстояние между концами отрезков */
			return std::min<T>(	{	
									norm(s1.point(0.0) - s2.point(0.0)), norm(s1.point(0.0) - s2.point(1.0)),
									norm(s1.point(1.0) - s2.point(0.0)), norm(s1.point(1.0) - s2.point(1.0)) 
								});
		}
		else
		{
			/* Искомое расстояние */
			return norm(s1.point(sc) - s2.point(tc));
		}
	}
}

/* Кейзы */
//#define CASE_OVERLAPPED
#ifndef CASE_OVERLAPPED
#define CASE_NOT_OVERLAPPED
#endif

//#define CASE_INTERSECTED
//#define CASE_PARALLEL
//#define CASE_COMPLANAR
//#define CASE_SKEW
//#define CASE_COLLINEAR
#define CASE_PERPENDICULAR

#if defined(CASE_INTERSECTED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 2.0, 0.0, 0.0 };

	/* Точки второго отрезка */
	auto Q0 = Point{ 1.0, -1.0, 1.0 };
	auto Q1 = Point{ 1.0, 1.0, -1.0 };

#elif defined(CASE_PARALLEL) && defined(CASE_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 6.0, 0.0, 0.0 };

	/* Точки второго отрезка */
	auto Q0 = Point{ 4.0, 5.0, 0.0 };
	auto Q1 = Point{ 8.0, 5.0, 0.0 };

#elif defined(CASE_PARALLEL) && defined(CASE_NOT_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 4.0, 0.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 5.0, 5.0, 0.0 };
	auto Q1 = Point{ 9.0, 5.0, 0.0 };

#elif defined(CASE_COMPLANAR) && defined(CASE_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 5.0, 0.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 7.0, 1.0, 0.0 };
	auto Q1 = Point{ 6.0, 4.0, 0.0 };

#elif defined(CASE_COMPLANAR) && defined(CASE_NOT_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ -5.0, 0.0, 0.0 };
	auto P1 = Point{ 0.0, 0.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 1.0, 2.0, 0.0 };
	auto Q1 = Point{ 5.0, 5.0, 0.0 };

#elif defined(CASE_SKEW) && defined(CASE_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 5.0, 5.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 5.0, 0.0, 1.0 };
	auto Q1 = Point{ 0.0, 5.0, 3.0 };

#elif defined(CASE_SKEW) && defined(CASE_NOT_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 5.0, 5.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 2.0, 3.0, 1.0 };
	auto Q1 = Point{ 0.0, 5.0, 3.0 };

#elif defined(CASE_COLLINEAR) && defined(CASE_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 5.0, 0.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 2.0, 0.0, 0.0 };
	auto Q1 = Point{ 8.0, 0.0, 0.0 };

#elif defined(CASE_COLLINEAR) && defined(CASE_NOT_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 5.0, 0.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 7.0, 0.0, 0.0 };
	auto Q1 = Point{ 6.0, 0.0, 0.0 };

#elif defined(CASE_PERPENDICULAR) && defined(CASE_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 5.0, 0.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 3.0, 3.0, 0.0 };
	auto Q1 = Point{ 3.0, 8.0, 0.0 };

#elif defined(CASE_PERPENDICULAR) && defined(CASE_NOT_OVERLAPPED)
	/* Точки первого отрезка */
	auto P0 = Point{ 0.0, 0.0, 0.0 };
	auto P1 = Point{ 5.0, 0.0, 0.0 };
	
	/* Точки второго отрезка */
	auto Q0 = Point{ 6.0, 1.0, 0.0 };
	auto Q1 = Point{ 6.0, 6.0, 0.0 };

#endif

int main()
{
	/* Задано два отрезка в пространстве */
	auto S1 = Segment{ P0, P1 };
	auto S2 = Segment{ Q0, Q1 };

	/* Необходимо найти расстояние между ними */
	print("Distance equals:", distance(S1, S2));
}
