#ifndef PLANT_EXTRINSIC_DRIVERS_H
#define PLANT_EXTRINSIC_DRIVERS_H

#include <plant/util.h>
#include <type_traits>

namespace plant {

//template<typename T, typename U>
//class PairVariant {
//		static_assert(!std::is_same<T, U>::value);
//public:
//		//PairVariant() = default;
//		PairVariant(T const& t_) {
//				data_.t = t_;
//				using_t = true;
//		}
//		PairVariant(U const& u_) {
//				data_.u = u_;
//				using_t = false;
//		}
//
//		template<typename K>
//		bool holds_alternative() { // ambiguous case? or do we prevent variants with the same type being constructed
//				if (std::is_same<K, T>::value && using_t) {
//						return true;
//				} else if (std::is_same<K, U>::value && !using_t) {
//						return true;
//				}
//				return false;
//		}
//
//		template<typename K>
//		const K& get() { // compiler may warn "returns reference to temporary" but its guaranteed to have active lifetime
//				if (holds_alternative<K>()) {
//						return using_t ? data_.t : data_.u;
//				} else {
//						util::stop("Tried to access type in variant that doesn't actively exist in variant.");
//				}
//		}
//
//private:
//		union Data {
//				T t;
//				U u;
//
//				~Data() {}; // not sure if I need this line
//		};
//		Data data_;
//		bool using_t;
//};

class Function {
public:
		Function(bool is_const);
private:
		interpolator::Interpolator variable;
		double constant;
		bool is_constant;
};

class ExtrinsicDrivers {
		using variable = interpolator::Interpolator;
public:
		// this will override any previously defined drivers with the same name
		void set_constant(std::string driver_name, double k) {
				PairVariant<double, interpolator::Interpolator> x = k;
				drivers[driver_name] = x;
		}

		// initialise spline of driver with x, y control points
		void set_variable(std::string driver_name, std::vector<double> const& x, std::vector<double> const& y) {
				drivers[driver_name] = variable();
				drivers[driver_name].get<variable>().init(x, y);
				set_extrapolation(driver_name, false); // default no extrapolation
		}

		void set_extrapolation(std::string driver_name, bool extrapolate) {
				drivers.at(driver_name).set_extrapolate(extrapolate);
		}

		// evaluate/query interpolated spline for driver at point u, return s(u), where s is interpolated function
		double evaluate(std::string driver_name, double u) const {
				return drivers.at(driver_name).eval(u);
		}

		// evaluate/query interpolated spline for driver at vector of points, return vector of values
		std::vector<double> evaluate_range(std::string driver_name, std::vector<double> u) const {
				return drivers.at(driver_name).r_eval(u);
		}

		// returns the name of each active driver - useful for R output
		std::vector<std::string> get_names() {
				auto ret = std::vector<std::string>();
				for (auto const& driver : drivers) {
						ret.push_back(driver.first);
				}
				return ret;
		}
private:
		std::unordered_map<std::string, Function> drivers;
};

}

#endif //PLANT_EXTRINSIC_DRIVERS_H
