both = %w(sin cos tan asin acos atan sinh cosh tanh asinh acosh atanh exp log sqrt)
fonly = %w(exp2 expm1 log10 log1p log2 logb cbrt erf erfc lgamma tgamma ceil floor round trunc sign)
dz = %w(real imag arg)

File.open('elementwise_function.c', 'w') do |f|
	f.puts 'static double sign(double x){ return ( x==0.0 ? 0.0 : ( x>0.0 ? 1.0 : ( x<0.0 ? -1.0 : x ))); }'
	
	both.each do |func|
		f.puts <<"EOS"
static void
km_#{func}_d(double *ent, void *null)
{
	*ent = #{func}(*ent);
}
static void
km_#{func}_z(COMPLEX *ent, void *null)
{
	*ent = c#{func}(*ent);
}
static void
km_#{func}_v(VALUE *ent, void *null)
{
	*ent = rb_funcall(rb_mMath, id_#{func}, 1, *ent);
}
VALUE
kmm_mat_#{func}_dest(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_DOUBLE ) {
		km_smat_each_d(smat, km_#{func}_d, NULL);
	} else if ( smat->vtype == VT_COMPLEX ) {
		km_smat_each_z(smat, km_#{func}_z, NULL);
	} else if ( smat->vtype == VT_VALUE ) {
		km_smat_each_v(smat, km_#{func}_v, NULL);
	} else {
		rb_raise(km_eVT, "the method is available only for float or complex matricies");
	}
	return self;
}
EOS
	end
	
	fonly.each do |func|
		f.puts <<"EOS"
static void
km_#{func}_d(double *ent, void *null)
{
	*ent = #{func}(*ent);
}
static void
km_#{func}_v(VALUE *ent, void *null)
{
	*ent = rb_funcall(rb_mMath, id_#{func}, 1, *ent);
}
VALUE
kmm_mat_#{func}_dest(VALUE self)
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype == VT_DOUBLE ) {
		km_smat_each_d(smat, km_#{func}_d, NULL);
	} else if ( smat->vtype == VT_VALUE ) {
		km_smat_each_v(smat, km_#{func}_v, NULL);
	} else {
		rb_raise(km_eVT, "the method is available only for float matricies");
	}
	return self;
}
EOS
	end
	f.puts <<'EOS'
static VALUE
km_zmat_e_funcapp(VALUE self, void (*func)(COMPLEX *,  void *))
{
	km_check_frozen(self);
	SMAT *smat = km_mat2smat(self);
	if ( smat->vtype != VT_COMPLEX ) {
		rb_raise(km_eVT, "the method is available only for complex matrcies");
	}
	km_smat_each_z(smat, func, NULL);
	return self;
}
static VALUE
km_zdmat_e_funcapp(VALUE self, VALUE op, void (*func)(double *, const COMPLEX *,  void *))
{
	km_check_frozen(op);
	SMAT *dest = km_mat2smat(op), *src = km_mat2smat(self);
	if ( dest->vtype != VT_DOUBLE || src->vtype != VT_COMPLEX ) {
		rb_raise(km_eVT, "self(operand) and argument(output) must be complex and float matrix, respectively");
	}
	CHECK_SAME_SIZE(dest, src);
	km_smat_each2_dcz(dest, src, func, NULL);
	return op;
}
EOS
	dz.each do |func|
		f.puts <<"EOS"
static void
km_#{func}_z(COMPLEX *ent, void *null)
{
	*ent = cpack(c#{func}(*ent), 0.0);
}
static void
km_#{func}_dcz(double *dest, const COMPLEX *src, void *null)
{
	*dest = c#{func}(*src);
}
VALUE
kmm_mat_#{func}_destl(int argc, VALUE *argv, VALUE self)
{
	rb_check_arity(argc, 0, 1);
	if ( argc == 0 ) {
		return km_zmat_e_funcapp(self, km_#{func}_z);
	} else {
		return km_zdmat_e_funcapp(self, argv[0], km_#{func}_dcz);
	}
}
VALUE
kmm_mat_#{func}(VALUE self)
{
	SMAT *smat = km_mat2smat(self);
	return km_zdmat_e_funcapp(self, km_Mat(smat->m, smat->n, VT_DOUBLE), km_#{func}_dcz);
}
EOS
	end
end

File.open('elementwise_function.h', 'w') do |f|
	(both+fonly).each do |func|
		f.puts "VALUE kmm_mat_#{func}_dest(VALUE self);"
	end
	dz.each do |func|
		f.puts "VALUE kmm_mat_#{func}_destl(int argc, VALUE *argv, VALUE self);"
		f.puts "VALUE kmm_mat_#{func}(VALUE self);"
	end
end
File.open('elementwise_function_definitions.c', 'w') do |f|
	f.puts "static void\nkm_define_efs(void)\n{"
	(both+fonly).each do |func|
		f.puts %Q|\trb_define_method(km_cMat, "#{func}!", kmm_mat_#{func}_dest, 0);|
		f.puts %Q|\trb_funcall(km_cMat, id__define_dup_escaped_method, 1, rb_str_new_cstr("#{func}"));|
	end
	dz.each do |func|
		f.puts %Q|\trb_define_method(km_cMat, "#{func}!", kmm_mat_#{func}_destl, -1);|
		f.puts %Q|\trb_define_method(km_cMat, "#{func}", kmm_mat_#{func}, 0);|
	end
	f.puts "}"
end
