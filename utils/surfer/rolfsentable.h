//
//  rolfsentable.h
//
//  GENERATED FILE -- do not edit by hand.
//  Produced by tests/gen_knot_census_names.cpp + tests/gen_knot_names_header.py
//  from tests/pd_codes_up_to_13_crossings.csv (KnotInfo), cross-validated against SnapPy.
//

#ifndef ROLFSENTABLE_H
#define ROLFSENTABLE_H

#include <optional>
#include <string>
#include <unordered_map>

/*! \file utils/surfer/rolfsentable.h
 *  \brief Translates Regina census hit names (e.g. "m004 : #1", "L104001")
 *  into classical Rolfsen knot table names (e.g. "4_1"), for knots up to
 *  10 crossings -- the largest crossing number Rolfsen's table covers.
 */

namespace rolfsen {

/**
 * Base census manifold name (the " : #N" suffix stripped) -> classical
 * Rolfsen name. Keyed by base name, not the full raw hit string: entries
 * sharing a base name in Regina's cusped hyperbolic census are different
 * ideal triangulations of the *same* manifold (verified by identical
 * volume across every such group), not different manifolds, so matching
 * on the base name is what makes this table robust to which specific
 * triangulation a caller's simplify() happens to land on.
 */
inline const std::unordered_map<std::string, std::string> table = {
    {"L104001", "4_1"},
    {"L105002", "5_2"},
    {"L106001", "6_1"},
    {"L107002", "7_2"},
    {"L107003", "7_3"},
    {"L107005", "7_5"},
    {"L108004", "8_4"},
    {"L108009", "8_9"},
    {"L108010", "8_10"},
    {"L108013", "8_13"},
    {"L108018", "8_18"},
    {"L108021", "8_21"},
    {"L109002", "9_2"},
    {"L109003", "9_3"},
    {"L109004", "9_4"},
    {"L109005", "9_5"},
    {"L109009", "9_9"},
    {"L109010", "9_10"},
    {"L109011", "9_11"},
    {"L109012", "9_12"},
    {"L109019", "9_19"},
    {"L109025", "9_25"},
    {"L109026", "9_26"},
    {"L109027", "9_27"},
    {"L109034", "9_34"},
    {"L109035", "9_35"},
    {"L109042", "9_42"},
    {"L109044", "9_44"},
    {"L109045", "9_45"},
    {"L109046", "9_46"},
    {"L109048", "9_48"},
    {"L109049", "9_49"},
    {"L110002", "10_2"},
    {"L110004", "10_4"},
    {"L110006", "10_6"},
    {"L110008", "10_8"},
    {"L110009", "10_9"},
    {"L110016", "10_16"},
    {"L110017", "10_17"},
    {"L110036", "10_36"},
    {"L110042", "10_42"},
    {"L110047", "10_47"},
    {"L110048", "10_48"},
    {"L110060", "10_60"},
    {"L110061", "10_61"},
    {"L110069", "10_69"},
    {"L110076", "10_76"},
    {"L110083", "10_86"},
    {"L110094", "10_94"},
    {"L110123", "10_123"},
    {"L110125", "10_125"},
    {"L110127", "10_127"},
    {"L110128", "10_128"},
    {"L110131", "10_131"},
    {"L110132", "10_132"},
    {"L110133", "10_133"},
    {"L110134", "10_134"},
    {"L110136", "10_136"},
    {"L110137", "10_137"},
    {"L110139", "10_139"},
    {"L110144", "10_144"},
    {"L110146", "10_146"},
    {"L110152", "10_152"},
    {"L110158", "10_158"},
    {"L110164", "10_163"},
    {"L110165", "10_164"},
    {"m004", "4_1"},
    {"m015", "5_2"},
    {"m032", "6_1"},
    {"m053", "7_2"},
    {"m074", "8_1"},
    {"m094", "9_2"},
    {"m199", "9_42"},
    {"m201", "10_132"},
    {"m222", "8_20"},
    {"m289", "6_2"},
    {"m340", "7_3"},
    {"m372", "9_46"},
    {"m389", "10_139"},
    {"o9_32208", "10_5"},
    {"o9_37080", "10_136"},
    {"o9_37732", "10_133"},
    {"o9_37770", "8_8"},
    {"o9_39277", "10_141"},
    {"o9_39339", "9_35"},
    {"o9_40064", "9_7"},
    {"o9_40076", "9_9"},
    {"o9_41611", "9_8"},
    {"o9_42258", "8_11"},
    {"o9_42277", "9_11"},
    {"o9_42320", "10_9"},
    {"o9_42963", "10_6"},
    {"o9_42974", "10_134"},
    {"o9_43592", "8_13"},
    {"o9_43609", "10_152"},
    {"o9_43771", "9_45"},
    {"o9_43874", "8_10"},
    {"o9_44057", "9_10"},
    {"s016", "10_1"},
    {"s385", "10_125"},
    {"s526", "8_2"},
    {"s558", "9_3"},
    {"s580", "10_145"},
    {"s648", "7_4"},
    {"s704", "10_140"},
    {"s726", "8_3"},
    {"s862", "8_4"},
    {"s870", "9_4"},
    {"s912", "6_3"},
    {"t09859", "10_142"},
    {"t09901", "10_130"},
    {"t10499", "10_126"},
    {"t10932", "8_5"},
    {"t11034", "8_7"},
    {"t11291", "7_6"},
    {"t11675", "9_6"},
    {"t12200", "10_153"},
    {"t12271", "9_44"},
    {"t12395", "8_6"},
    {"t12587", "8_9"},
    {"t12656", "7_7"},
    {"v1217", "10_2"},
    {"v2166", "10_161"},
    {"v2284", "9_5"},
    {"v2362", "10_3"},
    {"v2488", "10_4"},
    {"v2553", "10_128"},
    {"v2858", "10_8"},
    {"v3310", "7_5"},
    {"v3505", "8_21"},
};

/**
 * Returns the Rolfsen name (e.g. "4_1") corresponding to a raw Regina
 * census hit name (e.g. "m004 : #1" or "L104001"), or nullopt if unknown.
 * Any " : #N" suffix is stripped before looking up in table (see its
 * comment above), so callers can pass CensusHit::name() unmodified.
 */
inline std::optional<std::string> rolfsenName(const std::string &censusName) {
    std::string base = censusName.substr(0, censusName.find(" : "));
    auto it = table.find(base);
    if (it == table.end())
        return std::nullopt;
    return it->second;
}

} // namespace rolfsentable

#endif // ROLFSENTABLE_H 
