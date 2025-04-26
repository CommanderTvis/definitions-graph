const Area = {
    SetTheory: "Set Theory",
    OrderTheory: "Order Theory",
    GroupTheory: "Group Theory",
    RingTheory: "Ring Theory",
    FieldTheory: "Field Theory",
    NumberSystems: "Number Systems",
    Analysis: "Real Analysis",
    RealFunctions: "Real Functions",
    LinearAlgebra: "Linear Algebra",
    Geometry: "Geometry",
};

class Definition {
    constructor(area, id, name, refNames = [], url = null, importance = 0) {
        this.area = area;
        this.id = id;
        this.name = name;
        this.url = url;
        this.refNames = refNames;
        this.importance = importance;
    }

    get references() {
        return this.refNames.map(it => allDefinitions.filter(j => j.id === it)[0]);
    }
}

function proofwiki(name) {
    return "https://proofwiki.org/wiki/Definition:" + name;
}

const allDefinitions = [
    new Definition(Area.SetTheory, "set", "Set {}", [], proofwiki("Set"), 2),
    new Definition(Area.SetTheory, "subset", "Subset", ["set"], proofwiki("Subset")),
    new Definition(
        Area.SetTheory, "empty_set", "Empty Set ∅", ["set"], proofwiki("Empty_Set")),
    new Definition(Area.SetTheory, "union", "Union ∪", ["set"], proofwiki("Set_Union")),
    new Definition(Area.SetTheory, "power_set", "Power Set 𝒫 (A)", ["set"], proofwiki("Power_Set")),
    new Definition(
        Area.SetTheory, "intersection", "Intersection ∩", ["union"], proofwiki("Set_Intersection")),
    new Definition(
        Area.SetTheory,
        "cartesian_product",
        "Cartesian Product ×", ["tuple", "power_set", "union"], proofwiki("Cartesian_Product")),
    new Definition(
        Area.SetTheory, "tuple", "Ordered Pair (a, b)", ["set"], proofwiki("Ordered_Pair")),
    new Definition(
        Area.SetTheory,
        "disjoint_union",
        "Disjoint Union ⊔", ["union", "intersection"], proofwiki("Disjoint_Union")),
    new Definition(
        Area.SetTheory, "complement", "Set Difference \\\\", ["set"], proofwiki("Set_Difference")),
    new Definition(
        Area.SetTheory, "inductive_set", "Inductive Set", ["set"], proofwiki("Inductive_Set")),
    new Definition(
        Area.SetTheory,
        "function", "Mapping A → B", ["cartesian_product"], proofwiki("Mapping"), 2),
    new Definition(
        Area.SetTheory, "restriction", "Restriction (Mapping)", ["function"], proofwiki("Restriction/Mapping")),
    new Definition(Area.SetTheory, "injection", "Injection", ["function"], proofwiki("Injection/Definition_1")),
    new Definition(Area.SetTheory, "surjection", "Surjection", ["function"], proofwiki("Surjection/Definition_1")),
    new Definition(
        Area.SetTheory,
        "bijection", "Bijection", ["injection", "surjection"], proofwiki("Bijection/Definition_1"), 1),
    new Definition(Area.SetTheory, "image", "Image (Mapping)", ["function"], proofwiki("Image_(Set_Theory)/Mapping/Mapping")),
    new Definition(Area.SetTheory, "preimage", "Preimage (Mapping)", ["function"], proofwiki("Preimage/Mapping/Mapping")),
    new Definition(
        Area.SetTheory, "composition", "Composition ∘", ["function"], proofwiki("Composition_of_Mappings")),
    new Definition(
        Area.SetTheory, "permutation", "Permutation", ["bijection"], proofwiki("Permutation#:~:text=Definition-,A%20bijection,.,-Permutation%20on")),
    new Definition(
        Area.SetTheory,
        "sequence",
        "Sequence ℕ → A",
        ["function", "natural_numbers"], proofwiki("Sequence"), 1),
    new Definition(
        Area.SetTheory,
        "subsequence",
        "Subsequence ℕ → ℕ → A", ["sequence", "composition"], proofwiki("Subsequence")),
    new Definition(
        Area.SetTheory,
        "equivalence_relation",
        "Equivalence Relation", ["cartesian_product"], proofwiki("Equivalence_Relation/Definition_1"), 1),
    new Definition(
        Area.SetTheory,
        "equivalence_class",
        "Equivalence Class", ["equivalence_relation"], proofwiki("Equivalence_Class"), 1),
    new Definition(
        Area.SetTheory,
        "quotient_set",
        "Quotient Set", ["equivalence_class", "union", "intersection"], proofwiki("Quotient_Set")),
    new Definition(
        Area.SetTheory,
        "set_equivalence", "Set Equivalence", ["bijection"], proofwiki("Set_Equivalence")),
    new Definition(
        Area.SetTheory,
        "cardinal", "Cardinal", ["set_equivalence"], proofwiki("Cardinal")),

    new Definition(Area.OrderTheory,
        "order", "Ordering", ["cartesian_product"], proofwiki("Ordering"), 2),
    new Definition(Area.OrderTheory,
        "linear_order", "Total Ordering", ["order"], proofwiki("Total_Ordering"), 2),
    new Definition(
        Area.OrderTheory,
        "greatest_least_element", "Greatest Element", ["order"], proofwiki("Greatest_Element")),
    new Definition(
        Area.OrderTheory,
        "maximal_minimal_element", "Maximal Element", ["order"], proofwiki("Maximal_Element")),
    new Definition(Area.OrderTheory,
        "upper_lower_bound",
        "Upper Bound (Set)", ["order"], proofwiki("Upper_Bound_of_Set"), 1),
    new Definition(
        Area.OrderTheory,
        "supremum_infimum", "Supremum (Set)", ["upper_lower_bound"], proofwiki("Supremum of Set")),
    new Definition(
        Area.OrderTheory,
        "monotonic_function", "Monotone Mapping", ["order", "function"], proofwiki("Monotone_(Order_Theory)/Mapping")),
    new Definition(
        Area.OrderTheory,
        "monotonic_sequence",
        "Monotone Sequence", ["order", "sequence"], proofwiki("Monotone_Sequence")),
    new Definition(
        Area.OrderTheory,
        "bounded_function",
        "Bounded Mapping", ["function", "bounded_set", "image"], proofwiki("Bounded_Mapping")),
    new Definition(
        Area.OrderTheory,
        "bounded_set", "Bounded Set", ["order", "set"], proofwiki("Bounded_Set")),

    new Definition(
        Area.GroupTheory, "monoid", "Monoid", ["function"], proofwiki("Monoid"), 1),
    new Definition(
        Area.GroupTheory, "group", "Group", ["monoid"], proofwiki("Group"), 2),
    new Definition(
        Area.GroupTheory, "abelian_group", "Abelian Group", ["group"], proofwiki("Abelian_Group")),
    new Definition(
        Area.GroupTheory,
        "direct_group_product", "Group Direct Product", ["group", "cartesian_product"], proofwiki("Group Direct Product")),
    new Definition(
        Area.GroupTheory,
        "group_normalizer", "Group Normalizer", ["subgroup"], proofwiki("Normalizer")),
    new Definition(
        Area.GroupTheory,
        "group_centralizer", "Centralizer", ["group"], proofwiki("Centralizer")),
    new Definition(
        Area.GroupTheory,
        "group_center", "Group Center", ["group_centralizer"], proofwiki("Center_of_Group")),
    new Definition(
        Area.GroupTheory, "group_commutator", "Group Commutator", ["group"], proofwiki("Commutator")),
    new Definition(
        Area.GroupTheory, "subgroup", "Subgroup", ["group"], proofwiki("Subgroup"), 2),
    new Definition(
        Area.GroupTheory,
        "trivial_subgroup", "Trivial Subgroup", ["group"], proofwiki("Trivial_Subgroup")),
    new Definition(
        Area.GroupTheory,
        "conjugate_class", "Coset", ["subgroup"], proofwiki("Coset"), 1),
    new Definition(
        Area.GroupTheory,
        "quotient_group",
        "Quotient Group", ["conjugate_class"], proofwiki("Quotient_Group"), 1),
    new Definition(
        Area.GroupTheory,
        "subgroup_index",
        "Subgroup Index", ["subgroup", "conjugate_class"], proofwiki("Index_of_Subgroup")),
    new Definition(
        Area.GroupTheory, "group_commutant", "Derived Subgroup", ["subgroup", "group_commutator"], proofwiki("Derived_Subgroup")),
    new Definition(
        Area.GroupTheory,
        "normal_subgroup",
        "Normal Subgroup",
        ["subgroup"],
        proofwiki("Normal_Subgroup"),
        1),
    new Definition(Area.GroupTheory, "element_power", "Element Power", ["group"], proofwiki("Power_of_Element/Group")),
    new Definition(
        Area.GroupTheory, "cyclic_group", "Cyclic Group", ["element_power"], proofwiki("Cyclic_Group")),
    new Definition(
        Area.GroupTheory,
        "altering_group", "Alternating Group", ["symmetric_group", "homomorphism_kernel"], proofwiki("Alternating_Group")),
    new Definition(
        Area.GroupTheory,
        "symmetric_group",
        "Symmetric Group",
        ["group", "permutation", "composition"],
        proofwiki("Symmetric_Group")),
    new Definition(
        Area.GroupTheory,
        "group_homomorphism",
        "Group Homomorphism", ["function", "group"], proofwiki("Group_Homomorphism"), 1),
    new Definition(
        Area.GroupTheory,
        "group_epimorphism",
        "Group Epimorphism", ["group_homomorphism", "surjection"], proofwiki("Group_Epimorphism")),
    new Definition(
        Area.GroupTheory,
        "group_monomorphism",
        "Group Monomorphism", ["group_homomorphism", "injection"], proofwiki("Group_Monomorphism")),
    new Definition(
        Area.GroupTheory,
        "group_isomorphism",
        "Group Isomorphism", ["group_homomorphism", "bijection"], proofwiki("Group_Isomorphism"), 1),
    new Definition(
        Area.GroupTheory,
        "group_endomorphism", "Group Endomorphism", ["group_homomorphism"], proofwiki("Group_Endomorphism")),
    new Definition(
        Area.GroupTheory,
        "group_automorphism",
        "Group Automorphism", ["group_isomorphism"], proofwiki("Group_Automorphism")),
    new Definition(
        Area.GroupTheory,
        "homomorphism_kernel", "Homomorphism Kernel", ["group_homomorphism"], proofwiki("Kernel_of_Group_Homomorphism")),

    new Definition(
        Area.RingTheory, "ring", "Ring", ["abelian_group"], proofwiki("Ring_(Abstract_Algebra)"), 2),
    new Definition(
        Area.RingTheory,
        "commutative_ring", "Commutative Ring", ["ring"], proofwiki("Commutative_Ring")),
    new Definition(
        Area.RingTheory, "division_ring", "Division Ring", ["ring"], proofwiki("Division_Ring")),
    new Definition(
        Area.RingTheory,
        "integral_domain",
        "Integral Domain", ["commutative_ring", "zero_divisor"], proofwiki("Integral_Domain")),
    new Definition(
        Area.RingTheory, "zero_divisor", "Zero Divisor", ["ring"], proofwiki("Zero_Divisor/Ring")),
    new Definition(Area.RingTheory, "subring", "Subring", ["ring"], proofwiki("Subring")),
    new Definition(
        Area.RingTheory,
        "ideal", "Ring Ideal", ["subgroup", "ring"], proofwiki("Ideal_of_Ring"), 1),
    new Definition(
        Area.RingTheory, "maximal_ideal", "Maximal Ideal", ["ideal"], proofwiki("Maximal_Ideal_of_Ring")),
    new Definition(
        Area.RingTheory,
        "quotient_ring",
        "Quotient Ring", ["ring", "ideal", "quotient_group"], proofwiki("Quotient_Ring")),
    new Definition(
        Area.RingTheory, "prime_ideal", "Prime Ideal", ["ideal"], proofwiki("Prime_Ideal_of_Ring")),
    new Definition(
        Area.RingTheory, "principal_ideal", "Principal Ideal", ["ideal"], proofwiki("Principal_Ideal_of_Ring")),
    new Definition(
        Area.RingTheory,
        "principal_ideal_ring",
        "Principal Ideal Ring", ["principal_ideal", "ring"], proofwiki("Principal_Ideal_Ring")),
    new Definition(
        Area.RingTheory,
        "irreducible_element",
        "Irreducible Element", ["integral_domain"], proofwiki("Irreducible_Element_of_Ring")),
    new Definition(
        Area.RingTheory,
        "polynomial",
        "Polynomial", ["tuple", "commutative_ring"], proofwiki("Polynomial_over_Ring/One_Variable")),
    new Definition(
        Area.RingTheory,
        "polynomial_ring",
        "Polynomial Ring", ["polynomial", "commutative_ring"], proofwiki("Polynomial_Ring")),

    new Definition(Area.FieldTheory, "subfield", "Subfield", ["field"], proofwiki("Subfield")),
    new Definition(
        Area.FieldTheory,
        "field", "Field",
        ["division_ring"], proofwiki("Field_(Abstract_Algebra)/Definition_4"), 2),
    new Definition(Area.FieldTheory,
        "field_characteristic", "Field Characteristic", ["field", "integer_numbers"], proofwiki("Characteristic_of_Field"), 1),

    new Definition(
        Area.NumberSystems,
        "natural_numbers",
        "Natural Numbers ℕ", ["intersection", "inductive_set"], proofwiki("Natural_Numbers/Inductive_Sets_in_Real_Numbers"), 2),
    new Definition(
        Area.NumberSystems,
        "integer_numbers",
        "Integer Numbers ℤ", ["natural_numbers", "quotient_set"], proofwiki("Integer"), 2),
    new Definition(
        Area.NumberSystems,
        "rational_numbers",
        "Rational Numbers ℚ", ["integer_numbers"], proofwiki("Rational_Number"), 2),
    new Definition(
        Area.NumberSystems,
        "real_numbers",
        "Real Numbers ℝ", ["rational_numbers"], proofwiki("Real_Number"), 2),
    new Definition(
        Area.NumberSystems,
        "complex_numbers",
        "Complex Numbers ℂ", ["real_numbers"], proofwiki("Complex_Number"), 2),

    new Definition(
        Area.Analysis, "absolute_value", "Absolute Value", ["real_numbers"], proofwiki("Absolute_Value")),
    new Definition(
        Area.Analysis,
        "real_interval",
        "Numeric Interval", ["real_numbers"], proofwiki("Interval"), 1),
    new Definition(
        Area.Analysis,
        "real_sequence",
        "Numeric Sequence ℕ → ℝ", ["real_numbers", "sequence"], proofwiki("Sequence"), 1),
    new Definition(
        Area.Analysis,
        "real_function",
        "Real Function", ["real_numbers", "function"], proofwiki("Real_Function"), 2),
    new Definition(Area.RealFunctions, "sgn", "sgn", ["real_function"], proofwiki("Sign_Function")),
    new Definition(
        Area.Analysis,
        "delta_neighborhood",
        "Neighborhood", ["real_interval"], proofwiki("Neighborhood_(Real_Analysis)/Epsilon"), 2),
    new Definition(
        Area.Analysis,
        "accumulation_point",
        "Limit Point",
        ["real_numbers"], proofwiki("Limit_Point/Real_Analysis"), 1),
    new Definition(
        Area.Analysis,
        "function_limit",
        "Function Limit",
        ["accumulation_point", "real_function"],
        proofwiki("Limit_of_Real_Function"), 2),
    new Definition(
        Area.Analysis,
        "local_bounded_function",
        "Locally Bounded Function",
        ["accumulation_point", "real_function"],
        proofwiki("Locally_Bounded_Function")),
    new Definition(
        Area.Analysis,
        "infinitesimal_function",
        "Infinitesimal Function",
        ["real_function", "function_limit"],
        proofwiki("Infinitesimal")),
    new Definition(
        Area.Analysis,
        "continuous_at_point_function",
        "Continuous Function at Point",
        ["real_function", "delta_neighborhood", "intersection"],
        proofwiki("Continuous_Real_Function/Point"), 2),
    new Definition(
        Area.Analysis,
        "differentiable_at_point_function",
        "Function Differentiable at Point",
        ["real_function", "delta_neighborhood"],
        proofwiki("Differentiable_Mapping/Real_Function/Point"), 2),
    new Definition(
        Area.Analysis,
        "differential",
        "Function Differential at Point",
        ["differentiable_at_point_function", "real_function"],
        proofwiki("Differential"), 1),
    new Definition(
        Area.Analysis,
        "function_derivative",
        "Function Derivative at Point",
        ["function_limit"],
        proofwiki("Derivative"), 1),
    new Definition(
        Area.Analysis,
        "k_order_derivative",
        "k-th Order Derivative",
        ["function_derivative"],
        proofwiki("Higher_Derivative"), 1),
    new Definition(
        Area.Analysis,
        "one_sided_function_derivative",
        "Left and Right Derivatives",
        ["one_sided_limit"],
        proofwiki("One-Sided_Derivative")),
    new Definition(
        Area.Analysis,
        "uniform_continuous_function",
        "Uniformly Continuous Function",
        ["real_function", "absolute_value"],
        proofwiki("Uniform_Continuity")),
    new Definition(
        Area.Analysis,
        "continuous_function",
        "Continuous Function",
        ["continuous_at_point_function"],
        proofwiki("Continuous_Function"), 1),
    new Definition(
        Area.Analysis,
        "one_sided_limit",
        "One-sided Limit",
        ["function_limit"],
        proofwiki("One-Sided_Limit")),
    new Definition(
        Area.Analysis,
        "fundamental_sequence",
        "Cauchy Sequence",
        ["real_sequence", "delta_neighborhood"], proofwiki("Cauchy_Sequence")),
    new Definition(
        Area.Analysis,
        "deleted_delta_neighborhood",
        "Deleted Neighborhood",
        ["delta_neighborhood", "complement"], proofwiki("Deleted_Neighborhood/Real_Analysis"), 2),
    new Definition(
        Area.Analysis,
        "sequence_limit",
        "Sequence Limit",
        ["delta_neighborhood", "real_sequence"], proofwiki("Limit_of_Sequence"), 1),
    new Definition(
        Area.Analysis,
        "infinitesimal_sequence",
        "Infinitesimal Sequence", ["sequence_limit"], proofwiki("Infinitesimal_Sequence")),
    new Definition(
        Area.Analysis,
        "subsequential_limit",
        "Subsequential Limit", ["sequence_limit", "subsequence"], proofwiki("Subsequential_Limit")),
    new Definition(
        Area.Analysis,
        "limit_superior_inferior",
        "Limit Superior (Inferior)", ["supremum_infimum", "subsequential_limit"], proofwiki("Limit_Superior")),
    new Definition(
        Area.Analysis,
        "discontinuity", "Discontinuity", ["continuous_at_point_function"], proofwiki("Discontinuity_(Real_Analysis)")),
    new Definition(
        Area.RealFunctions,
        "exponential_function",
        "Exponential Function", ["real_function"], proofwiki("Exponential_Function")),
    new Definition(
        Area.RealFunctions,
        "logarithmic_function",
        "Logarithmic Function", ["real_function"], proofwiki("Logarithm")),
    new Definition(
        Area.RealFunctions,
        "power_function", "Power Function", ["real_function"], proofwiki("Power_Function")),
    new Definition(
        Area.RealFunctions, "sin_cos", "Sine and Cosine", ["real_function"], proofwiki("Sine")),
    new Definition(
        Area.RealFunctions, "tan_cot", "Tangent and Cotangent", ["sin_cos"], proofwiki("Tangent")),
    new Definition(
        Area.RealFunctions, "asin_acos", "Arcsine and Arccosine", ["sin_cos"], proofwiki("Arcsine")),
    new Definition(
        Area.RealFunctions, "atan_acot", "Arctangent and Arccotangent", ["tan_cot"], proofwiki("Arctangent")),
    new Definition(
        Area.Analysis,
        "small_o",
        "Little-o Notation",
        ["real_function", "deleted_delta_neighborhood", "absolute_value", "intersection"],
        proofwiki("Little-O_Notation"), 1),
    new Definition(
        Area.Analysis,
        "big_o",
        "Big-O Notation",
        ["real_function", "deleted_delta_neighborhood", "absolute_value", "intersection"],
        proofwiki("Big-O_Notation"), 1),
    new Definition(
        Area.Analysis,
        "extremum",
        "Extremum Point",
        ["real_function", "deleted_delta_neighborhood"], proofwiki("Extremum"), 1),
    new Definition(
        Area.Analysis,
        "parametric_curve",
        "Parametric Plane Curve", ["real_function"], proofwiki("Parametric_Equation")),
    new Definition(
        Area.Analysis,
        "series",
        "Series", ["tuple", "real_sequence"], proofwiki("Series")),
    new Definition(
        Area.Analysis,
        "convergent_series",
        "Convergent Series", ["series", "sequence_limit"], proofwiki("Convergent_Series")),
    new Definition(
        Area.Analysis,
        "convex_set",
        "Convex Numeric Set", ["real_numbers", "real_interval"], proofwiki("Convex_Set")),
    new Definition(
        Area.Analysis,
        "convex_function",
        "Convex Function", ["convex_set", "real_function"], proofwiki("Convex_Function")),

    new Definition(
        Area.LinearAlgebra,
        "linear_space",
        "Vector Space", ["abelian_group", "field"], proofwiki("Vector_Space"), 2),
    new Definition(
        Area.LinearAlgebra,
        "linear_subspace",
        "Vector Subspace", ["linear_space"], proofwiki("Subspace"), 1),
    new Definition(
        Area.LinearAlgebra,
        "linear_combination",
        "Linear Combination", ["linear_space"], proofwiki("Linear_Combination"), 1),
    new Definition(
        Area.LinearAlgebra,
        "linear_independence",
        "Linear Independence",
        ["linear_combination"], proofwiki("Linear_Independence"), 1),
    new Definition(
        Area.LinearAlgebra,
        "linear_span",
        "Linear Span", ["linear_space", "linear_combination"], proofwiki("Linear_Span")),
    new Definition(Area.LinearAlgebra,
        "linear_manifold", "Linear Manifold", ["linear_subspace"], proofwiki("Linear_Manifold")),
    new Definition(
        Area.LinearAlgebra,
        "basis",
        "Basis", ["linear_combination", "linear_space"], proofwiki("Basis_of_Vector_Space"), 2),
    new Definition(
        Area.LinearAlgebra,
        "orthogonal_basis",
        "Orthogonal Basis", ["basis", "dot_product"], proofwiki("Orthogonal_Basis"), 1),
    new Definition(
        Area.LinearAlgebra,
        "orthonormal_basis",
        "Orthonormal Basis",
        ["orthogonal_basis", "dot_product"], proofwiki("Orthonormal_Basis"), 1),
    new Definition(
        Area.LinearAlgebra,
        "dot_product",
        "Dot Product", ["linear_space"], proofwiki("Scalar_Product"), 1),
    new Definition(
        Area.LinearAlgebra,
        "cross_product",
        "Cross Product", ["linear_space", "orthonormal_basis"], proofwiki("Vector_Product")),
    new Definition(
        Area.LinearAlgebra,
        "triple_product", "Triple Product", ["cross_product", "dot_product"], proofwiki("Triple_Product")),

    new Definition(
        Area.LinearAlgebra,
        "matrix",
        "Matrix (m × n) → X", ["tuple", "natural_numbers"], proofwiki("Matrix"), 2),
    new Definition(
        Area.LinearAlgebra, "row_column", "Row (Column)", ["matrix"], proofwiki("Row_(Matrix)")),
    new Definition(Area.LinearAlgebra,
        "scalar_matrix", "Scalar Matrix", ["matrix"], proofwiki("Scalar_Matrix")),
    new Definition(
        Area.LinearAlgebra,
        "real_matrix_space",
        "Matrix Space", ["matrix", "linear_space", "field"], proofwiki("Matrix_Space"), 1),
    new Definition(
        Area.LinearAlgebra,
        "matrix_multiplication",
        "Matrix Multiplication",
        ["matrix", "field"], proofwiki("Matrix_Multiplication"), 1),
    new Definition(
        Area.LinearAlgebra,
        "elementary_matrix", "Elementary Matrix", ["matrix_multiplication"], proofwiki("Elementary_Matrix")),
    new Definition(Area.LinearAlgebra, "inverse_matrix", "Inverse Matrix", ["matrix_multiplication"], proofwiki("Inverse_Matrix")),
    new Definition(
        Area.LinearAlgebra,
        "determinant", "Determinant", ["permutation", "matrix"], proofwiki("Determinant"), 1),
    new Definition(Area.LinearAlgebra, "minor", "Minor", ["determinant"], proofwiki("Minor")),
    new Definition(Area.LinearAlgebra,
        "algebraic_complement", "Cofactor", ["minor"], proofwiki("Cofactor")),
    new Definition(
        Area.LinearAlgebra, "adjugate_matrix", "Adjugate Matrix", ["algebraic_complement"], proofwiki("Adjugate_Matrix")),
    new Definition(
        Area.LinearAlgebra, "rank", "Rank", ["linear_independence"], proofwiki("Rank_of_Matrix"), 1),
    new Definition(Area.LinearAlgebra, "linear_system", "System of Linear Equations", ["matrix"], proofwiki("System_of_Linear_Equations")),
    new Definition(Area.LinearAlgebra, "block_matrix", "Block Matrix", ["matrix"], proofwiki("Block_Matrix")),
    new Definition(Area.LinearAlgebra, "matrix_trace", "Matrix Trace", ["matrix"], proofwiki("Trace_of_Matrix"))
];

(() => {
    let arr = allDefinitions.map(it => it.id);
    for (const item of arr) {
        if (arr.indexOf(item) !== arr.lastIndexOf(item)) {
            throw Error(`Duplicate definition id: ${item}.`);
        }
    }

    for (const def of allDefinitions) {
        for (const ref of def.refNames) {
            if (!allDefinitions.some(it => it.id === ref)) {
                throw Error(`Unresolved reference to ${ref} in ${def.id}.`);
            }
        }
    }
})();
