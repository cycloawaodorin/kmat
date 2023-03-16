class Mat	
	def ===(other)
		if other.kind_of?(Mat)
			self == other
		else
			false
		end
	end
	
	def [](*idx)
		if block_given?
			tmp = self.bracket(*idx)
			ret = yield(tmp)
			tmp.__send__(:_kill)
			ret
		else
			self.bracket(*idx)
		end
	end
	
	def temp
		ret = yield(self)
		_kill
		ret
	end
	
	# for example, A.repmat(2, 3) returns [A, A, A: A, A, A]
	def repmat(row_repeat, col_repeat)
		ri = Mat.new(row_size()*row_repeat, 1, :int) do |i, j|
			i % row_size()
		end
		ci = Mat.new(col_size()*col_repeat, 1, :int) do |i, j|
			i % col_size()
		end
		self[ri, ci] do |m|
			m.dup
		end
	end
	
	# extend `other' to fit the shape of self
	# this is similar to Python's numpy broadcasting
	def broadcast(other, vt=nil)
		if other.kind_of?(Mat)
			if other.row_size == 1
				if other.col_size == 1
					other.repmat(row_size(), col_size())
				elsif other.col_size == self.col_size
					other.repmat(self.row_size, 1)
				else
					raise MismatchedDimensionError, "can't broadcast from (#{other.size.join(', ')}) to (#{self.size.join(', ')})"
				end
			elsif other.row_size == self.row_size
				if other.col_size == 1
					other.repmat(1, col_size())
				elsif other.col_size == self.col_size
					other
				else
					raise MismatchedDimensionError, "can't broadcast from (#{self.size.join(', ')}) to (#{other.size.join(', ')})"
				end
			else
				raise MismatchedDimensionError, "can't broadcast from (#{self.size.join(', ')}) to (#{other.size.join(', ')})"
			end
		else
			vt = vtype() if vt.nil?
			Mat.new(row_size(), col_size(), vt).fill(other)
		end
	end
	
	# replace self by repeated `src'
	# self will be left-top of infinitly repeated `src'
	def replace_rep(src)
		ri = Mat.new(row_size(), 1, :int) do |i, j|
			i % src.row_size
		end
		ci = Mat.new(col_size(), 1, :int) do |i, j|
			i % src.col_size
		end
		src[ri, ci] do |m|
			self.copy_from(m)
		end
	end
	
	# flip row, column or both axis
	# `dim' Symbol specifies fliping axis (:row, :col, :both, :none are available)
	def flip(dim=:both)
		dim = dim.to_sym if dim.respond_to?(:to_sym)
		m, n = *shape()
		case dim
		when :both, :b
			ri = Mat.new(m, 1, :int) do |i, j|
				m-i-1
			end
			ci = Mat.new(n, 1, :int) do |i, j|
				n-i-1
			end
			self[ri, ci]
		when :row, :r
			ri = Mat.new(m, 1, :int) do |i, j|
				m-i-1
			end
			self[ri, nil]
		when :col, :c
			ci = Mat.new(n, 1, :int) do |i, j|
				n-i-1
			end
			self[nil, ci]
		when :none, :n
			self[]
		else
			raise ArgumentError, "unknown axis symbol #{dim.inspect}"
		end
	end
	
	def diag(k=0)
		_diag_ul(k)
	end
end

class << Mat
	# for example, Mat.blocks([a, b], [c, d]) returns a single matrix of [a, b; c, d]
	def blocks(args)
		vstack(*(args.map do |row|
			hstack(*row)
		end))
	end
	
	def vstack(*mats)
		m, n = mats[0].shape
		1.upto(mats.size-1) do |i|
			mat = mats[i]
			raise MismatchedDimensionError, 'column-sizes must be the same' if mat.col_size != n
			m += mat.row_size
		end
		ret = Mat.new(m, n, mats[0].vtype)
		k = 0;
		n = Mat.irange(n)
		mats.each do |mat|
			ret[Mat.new(mat.row_size, 1, :int) do |i, j|
				i+k
			end, n] = mat
			k += mat.row_size
		end
		ret
	end
	
	def hstack(*mats)
		m, n = mats[0].shape
		1.upto(mats.size-1) do |i|
			mat = mats[i]
			raise MismatchedDimensionError, 'row-sizes must be the same' if mat.row_size != m
			n += mat.col_size
		end
		ret = Mat.new(m, n, mats[0].vtype)
		k = 0
		m = Mat.irange(m)
		mats.each do |mat|
			ret[m, Mat.new(mat.col_size, 1, :int) do |i, j|
				i+k
			end] = mat
			k += mat.col_size
		end
		ret
	end
end
