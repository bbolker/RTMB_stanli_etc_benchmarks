// Stan's lower/upper-bound transform and its log-Jacobian, so TMB objectives
// can be written on Stan's unconstrained scale (lower = 0 is just exp(u),
// log-Jacobian u).
template <class Type>
Type lub_trans(Type u, Type a, Type b) { return a + (b - a) * invlogit(u); }

template <class Type>
Type lub_logjac(Type u, Type a, Type b) {
  return log(b - a) - logspace_add(Type(0), -u) - logspace_add(Type(0), u);
}
