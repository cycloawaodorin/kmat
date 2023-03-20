class Random
	private def _kmat_rand_d(arg=nil)
		if arg
			self.rand(arg).to_f
		else
			self.rand()
		end
	end
	private def _kmat_rand_z(arg=nil)
		if arg
			Complex(self.rand(arg), self.rand(arg))
		else
			Complex(self.rand(), self.rand())
		end
	end
	private def _kmat_rand_i(arg=nil)
		case arg
		when NilClass
			self.rand(2)
		when Integer, Range
			self.rand(arg)
		when Float
			if 0.0 <= arg && arg <= 1.0
				( self.rand() < arg ) ? 1 : 0
			else
				raise ArgumentError, 'Float argument must be in [0, 1] for vt:int'
			end
		else
			raise ArgumentError, 'for vt:int, arg must be Integer(Uniform(0...arg)), Range(Uniform(arg)), Float(1 with probability arg, 0 with probability 1-arg) or omitted'
		end
	end
	private def _kmat_rand_b(arg=nil)
		if arg.nil?
			if self.rand(2) == 0
				false
			else
				true
			end
		elsif arg.kind_of?(Numeric) && 0.0 <= arg && arg <= 1.0
			( self.rand() < arg )
		else
			raise ArgumentError, 'for vt:bool, arg must be Float(true with probability arg, false with probability 1-arg) or omitted'
		end
	end
	private def _kmat_rand_v(arg=nil)
		if arg
			self.rand(arg)
		else
			self.rand()
		end
	end
	def kmat_vt_rand(vt, arg=nil)
		case Mat.vt_refine(vt)
		when :float
			_kmat_rand_d(arg)
		when :complex
			_kmat_rand_z(arg)
		when :int
			_kmat_rand_i(arg)
		when :bool
			_kmat_rand_b(arg)
		when :object
			_kmat_rand_v(arg)
		else
			raise Mat::InternalError, 'unknown value type'
		end
	end
end

class << Random
	using(Module.new do
		refine Random.singleton_class do
			alias :kmat_srand_org :srand
		end
	end)
	def srand(*args)
		@kmat_stored = false
		kmat_srand_org(*args)
	end
	def randn(*args)
		$MatRandom.randn(*args)
	end
end

module Kernel
	module_function def randn(*args)
		Random.randn(*args)
	end
end

class Mat
	def rand(arg=nil, random: $MatRandom)
		_rand0(random, arg)
	end
	def randn(random: $MatRandom)
		_randn0(random)
	end
end
class << Mat
	def rand(m=1, n=1, random: $MatRandom)
		_rand0(m, n, random)
	end
	def randn(m=1, n=1, random: $MatRandom)
		_randn0(m, n, random)
	end
end
