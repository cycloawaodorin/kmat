class Mat
	def all?(*args)
		if block_given?
			each do |ent|
				return false unless yield(ent, *args)
			end
		elsif args.size == 0
			if self.vtype == :bool
				each do |ent|
					return false unless ent
				end
			else
				return enum_for(__method__)
			end
		else
			each do |ent|
				args.each do |arg|
					return false unless arg === ent
				end
			end
		end
		return true
	end
	def any?(*args)
		if block_given?
			each do |ent|
				return true if yield(ent, *args)
			end
		elsif args.size == 0
			if self.vtype == :bool
				each do |ent|
					return true if ent
				end
			else
				return enum_for(__method__)
			end
		else
			each do |ent|
				args.each do |arg|
					return true if arg === ent
				end
			end
		end
		return false
	end
	
	def ==(other)
		return true if self.equal?(other)
		return false unless other.kind_of?(Mat)
		if self.vtype == other.vtype && self.size == other.size
			each_with_index do |e, i, j|
				return false unless e == other[i, j]
			end
			true
		else
			false
		end
	end
	def <(other)
		return false if self.equal?(other)
		unless other.kind_of?(Mat)
			raise ArgumentError, "Mat < #{other.class} is not defined"
		end
		if self.vtype != other.vtype
			raise ValueTypeError, "value types must be the same"
		elsif self.size != other.size
			raise MismatchedDimensionError, "the sizes must be the same, (#{self.size.join(', ')}) != (#{other.size.join(', ')})"
		end
		each_with_index do |e, i, j|
			return false unless e < other[i, j]
		end
		true
	end
	def <=(other)
		return false if self.equal?(other)
		unless other.kind_of?(Mat)
			raise ArgumentError, "Mat <= #{other.class} is not defined"
		end
		if self.vtype != other.vtype
			raise ValueTypeError, "value types must be the same"
		elsif self.size != other.size
			raise MismatchedDimensionError, "the sizes must be the same, (#{self.size.join(', ')}) != (#{other.size.join(', ')})"
		end
		each_with_index do |e, i, j|
			return false unless e <= other[i, j]
		end
		true
	end
	def >(other)
		return false if self.equal?(other)
		unless other.kind_of?(Mat)
			raise ArgumentError, "Mat > #{other.class} is not defined"
		end
		if self.vtype != other.vtype
			raise ValueTypeError, "value types must be the same"
		elsif self.size != other.size
			raise MismatchedDimensionError, "the sizes must be the same, (#{self.size.join(', ')}) != (#{other.size.join(', ')})"
		end
		each_with_index do |e, i, j|
			return false unless e > other[i, j]
		end
		true
	end
	def >=(other)
		return false if self.equal?(other)
		unless other.kind_of?(Mat)
			raise ArgumentError, "Mat > #{other.class} is not defined"
		end
		if self.vtype != other.vtype
			raise ValueTypeError, "value types must be the same"
		elsif self.size != other.size
			raise MismatchedDimensionError, "the sizes must be the same, (#{self.size.join(', ')}) != (#{other.size.join(', ')})"
		end
		each_with_index do |e, i, j|
			return false unless e >= other[i, j]
		end
		true
	end
	
	def eq(other)
		unless other.kind_of?(Mat) && other.size == self.size
			other = broadcast(other)
		end
		Mat.new(row_size(), col_size(), :bool).eq!(self, other)
	end
	def lt(other)
		unless other.kind_of?(Mat) && other.size == self.size
			other = broadcast(other)
		end
		Mat.new(row_size(), col_size(), :bool).lt!(self, other)
	end
	def le(other)
		unless other.kind_of?(Mat) && other.size == self.size
			other = broadcast(other)
		end
		Mat.new(row_size(), col_size(), :bool).le!(self, other)
	end
	def gt(other)
		unless other.kind_of?(Mat) && other.size == self.size
			other = broadcast(other)
		end
		Mat.new(row_size(), col_size(), :bool).gt!(self, other)
	end
	def ge(other)
		unless other.kind_of?(Mat) && other.size == self.size
			other = broadcast(other)
		end
		Mat.new(row_size(), col_size(), :bool).ge!(self, other)
	end
end
