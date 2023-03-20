class Mat
	def least_squares!(a, b, tor: nil, weight: nil, transpose: false)
		if tor
			raise ArgumentError, "tor and weight must not specify at once" if weight
			raise ArgumentError, "for conjugate LS, auto-transpose is not available" if transpose
			_ls_conj(a, b, tor)
		elsif weight
			raise ArgumentError, "for weighted LS, auto-transpose is not available" if transpose
			wls!(a, b, weight)
		else
			if transpose
				tls!(a, b)
			else
				_ls(a, b)
			end
		end
	end
	def least_squares(b, tor: nil, weight: nil, transpose: false)
		if tor
			raise ArgumentError, "tor and weight must not specify at once" if weight
			raise ArgumentError, "for conjugate LS, auto-transpose is not available" if transpose
			ls_conj(b, tor)
		elsif weight
			raise ArgumentError, "for weighted LS, auto-transpose is not available" if transpose
			wls(b, weight)
		else
			if transpose
				tls(b)
			else
				Mat.new(*b.shape, :float).__send__(:_ls, self, b)
			end
		end
	end
	alias :ls :least_squares
	alias :ls! :least_squares!

	def ls_conj!(a, b, tor: Float::EPSILON)
		_ls_conj(a, b, tor)
	end
	def ls_conj(b, tor: Float::EPSILON)
		Mat.new(self.col_size, 1, :float).ls_conj!(self, b, tor: tor)
	end
	
	# weighted least squares
	# weight `w' is a vector consists of standard deviations of errors
	# a; (m, n)-matrix
	# b: m-vector
	# x: n-vector
	def wls!(a, b, w)
		nw = w.length
		diag = Mat.new(nw, nw, :f)
		diag.diag.temp do |dd|
			begin
				dd.copy_from(w)
			rescue MismatchedDimensionError
				dd.tcopy_from(w)
			end
			dd.sqrt!
		end
		glm!(a, diag, b)
	end
	def wls(b, w)
		Mat.new(self.col_size, 1, :f).wls!(self, b, w)
	end
	
	def rand_orth(random: $MatRandom)
		_rand_orth(random)
	end
	def self.rand_orth(n, random: $MatRandom)
		Mat.new(n, n, :float).rand_orth(random: random)
	end
	
	def under(other)
		if square?
			begin
				solve(other)
			rescue UncomputableMatrixError
				ls(other)
			end
		else
			ls(other)
		end
	end
	def under!(a, b)
		if a.square?
			begin
				solve!(a, b)
			rescue UncomputableMatrixError
				ls!(a, b)
			end
		else
			ls!(a, b)
		end
	end
	def over(other)
		if other.square?
			begin
				other.tsolve(self).t!
			rescue UncomputableMatrixError
				other.tls(self).t!
			end
		else
			other.tls(self).t!
		end
	end
	def over!(a, b)
		if b.square?
			begin
				tsolve!(b, a).t!
			rescue UncomputableMatrixError
				tls!(b, a).t!
			end
		else
			tls!(b, a).t!
		end
	end
	def pinv
		n = self.size.min
		i = Mat.I(n)
		if square?
			begin
				solve(i)
			rescue UncomputableMatrixError
				ls(i)
			end
		else
			ls(i)
		end
	end
	def eigen_values
		if symmetry?
			symmetrize().sym_eigen_values()
		else
			ge_eigen_values()
		end
	end
	def evd
		if symmetry?
			symmetrize().sym_evd()
		else
			ge_evd()
		end
	end
	
	# inner product
	# let x and y are column vectors and M is a square matrix
	# the following is available
	# (1) x.iprod => x'*x
	# (2) x.iprod(y) => x'*y
	# (3) M.iprod(x) => x'*M*x
	# (4) M.iprod(x, y) => x'*M*y
	# if `inverse' is true, M is replaced by its Moore–Penrose inverse
	# in case (1) or (2), `inverse' is ignored
	# if x or y is row vector, take transpose automatically
	def iprod(*args, inverse: false)
		case  args.size
		when 0
			if vector? # case 1
				_iprod(self)
			else
				raise MismatchedDimensionError, "self must be a vector"
			end
		when 1
			o = args[0]
			if o.vector?
				if vector? # case 2
					_iprod(o)
				else # case 3
					if !square?
						raise MismatchecDimensionError, "M must be square for M.iprod(x)"
					elsif inverse
						temp = Mat.new(row_size, 1, vtype)
						if o.col_size == 1
							begin
								temp.solve!(self, o)
							rescue UncomputableMatrixError
								temp.ls!(self, o)
							end
						else
							begin
								temp.tsolve!(self, o)
							rescue UncomputableMatrixError
								temp.tls!(self, o)
							end
						end
					else
						temp = Mat.new(row_size, 1, vtype)
						if o.col_size == 1
							temp.mprod!(self, o)
						else
							temp.mprod!(o, self)
						end
					end
					temp.__send__(:_iprod, o)
				end
			else
				raise MismatchedDimensionError, "as iprod with single argument, x.iprod(y) and M.iprod(x) are available, not x.iprod(M) or others"
			end
		when 2 # case 4
			x, y = *args
			if !square? || !x.vector? || !y.vector? || x.length != row_size || y.length != row_size
				raise MismatchecDimensionError, "as iprod with 2 arguments, only M.iprod(x, y) is available, not others"
			end
			temp = Mat(row_size, 1, vtype)
			if y.col_size == 1
				temp.mprod!(y)
			else
				if y.frozen?
					y = y.dup
					y.transpose
					temp.mprod!(y)
				else
					begin
						y.transpose
						temp.mprod!(y)
					ensure
						y.transpose
					end
				end
			end
			x.__send__(:_iprod, temp)
		else
			raise ArgumentError, "wrong number of arguments (given #{args.size}, expected 0..2)"
		end
	end
	alias :inner_product :iprod
	
	# distance (Frobenius norm of the difference)
	def distance(other)
		self.dup.sub!(other).normf
	end
	def squared_distance(other)
		self.dup.sub!(other).iprod
	end
	
	# matrix norm
	# `type' Symbol or Integer specifies norm definition
	# if `elementwise' is true, return induced norm
	# if `elementwise' is false, return elementwise norm
	# if `type' is :fro, return Frobenius norm
	def norm(type=:two, elementwise: false)
		if elementwise
			case type
			when :one, :o, 1
				norm_e1()
			when :infinity, :i, :inf, :infty, Float::INFINITY
				norm_einf()
			when :two, :t, 2, :frobenius, :f, :fro
				normf()
			end
		else
			case type
			when :one, :o, 1
				norm1()
			when :infinity, :i, :inf, :infty, Float::INFINITY
				normi()
			when :two, :t, 2
				norm2()
			when :frobenius, :f, :fro
				normf()
			end
		end
	end
	
	def rcond(type=:two)
		case type
		when :one, :o, 1
			_rcondoi(:one)
		when :infinity, :i, :inf, :infty, Float::INFINITY
			_rcondoi(:infinity)
		when :two, :t, 2
			_rcond2()
		when :frobenius, :f, :fro
			_rcondf()
		end
	end
end
