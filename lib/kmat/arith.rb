class Mat
	# to define Numeric * Mat as Mat#scalar(Numeric) remaining Mat * Mat is undefined,
	# this returns non‐standard value
	def coerce(other)
		if caller.first[/`([^']*)'/, 1] == '*'
			[self, other]
		else
			raise TypeError, "Mat can't be coerced into #{other.class}"
		end
	end
	def add(other)
		if other.kind_of?(Mat)
			self.dup.add!(broadcast(other))
		else
			self.dup.s_add!(other)
		end
	end
	alias :+ :add
	def sub(other)
		if other.kind_of?(Mat)
			self.dup.sub!(broadcast(other))
		else
			self.dup.s_sub!(other)
		end
	end
	alias :- :sub
	def e_mul(other)
		if other.kind_of?(Mat)
			self.dup.e_mul!(broadcast(other))
		else
			self.dup.s_mul!(other)
		end
	end
	def e_div(other)
		if other.kind_of?(Mat)
			self.dup.e_div!(broadcast(other))
		else
			self.dup.s_div!(other)
		end
	end
	module MatrixProductOperator
		refine Mat do
			def *(other)
				if other.kind_of?(Mat)
					mprod(other)
				else
					self.dup.s_mul!(other)
				end
			end
			def /(other)
				if other.kind_of?(Mat)
					over(other)
				else
					self.dup.s_div!(other)
				end
			end
		end
	end
	def *(other)
		if other.kind_of?(Mat)
			raise ArgumentError, "Mat#* is available only for scalar multiplication. To use Mat#*(Mat) for matrix product, call `using Mat::MatrixProductOperator'"
		else
			self.dup.s_mul!(other)
		end
	end
	def scalar(alpha)
		self.dup.s_mul!(alpha)
	end
	alias :"scalar!" :"s_mul!"
	def /(other)
		if other.kind_of?(Mat)
			raise ArgumentError, "Mat#/ is available only for scalar multiplication with reciprocal. To use Mat#/(Mat) as an alias of Mat#over, call `using Mat::MatrixProductOperator`"
		else
			self.dup.s_div!(other)
		end
	end
	
	def add_times(other, alpha)
		if other.kind_of?(Mat)
			self.dup.add_times!(broadcast(other), alpha)
		else
			self.dup.s_mul!(other*alpha)
		end
	end
	
	def maximum(other)
		if other.kind_of?(Mat)
			self.dup.maximum!(broadcast(other))
		else
			self.dup.s_maximum!(other)
		end
	end
	def minimum(other)
		if other.kind_of?(Mat)
			self.dup.minimum!(broadcast(other))
		else
			self.dup.s_minimum!(other)
		end
	end
	
	def pow(other)
		if other.kind_of?(Mat)
			self.dup.pow!(broadcast(other))
		else
			self.dup.s_pow!(other)
		end
	end
	def rpow(other)
		other.dup.pow!(other.broadcast(self))
	end
	def **(other)
		case other
		when Integer
			ret, temp = self.dup, self.dup
			ret.eye
			if 0 <= other
				a = self.dup
			else
				a = self.inv
				other = -other
			end
			loop do
				if other % 2 == 1
					temp.mprod!(ret, a)
					temp, ret = ret, temp
				end
				other = other.div(2)
				break if other == 0
				temp.mprod!(a, a)
				temp, a = a, temp
			end
			ret
		when Float
			if self.symmetry?
				v, d = self.symmetrize.sym_evd
				d.diag.s_pow!(other)
				foo = v.mprod(d)
				d.mprod!(foo, v.t!)
			else
				raise UncomputableMatrixError, "cannot compute float power of non-symmetric matrcies"
			end
		when Rational
			self ** other.to_f
		else
			raise ArgumentError, "Mat powered by #{other.class} is not supported"
		end
	end
	
	def hypot(other)
		if other.kind_of?(Mat)
			self.dup.hypot!(broadcast(other))
		else
			self.dup.s_hypot!(other)
		end
	end
end

class << Mat
	def max(*args)
		case args.size
		when 0
			nil
		when 1
			args[0]
		else
			ret = args[0].dup
			1.upto(args.size-1) { |i|
				ret.maximum!(args[i])
			}
			ret
		end
	end
	alias :maximum :max
	def min(*args)
		case args.size
		when 0
			nil
		when 1
			args[0]
		else
			ret = args[0].dup
			1.upto(args.size-1) { |i|
				ret.minimum!(args[i])
			}
			ret
		end
	end
	alias :minimum :min
end
