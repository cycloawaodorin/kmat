class Mat
	def sum(arg=:all)
		case arg
		when Mat
			_sum(arg)
		when :all, :a
			_sum(Mat.new(1, 1, self.vtype))
		when :col, :c
			_sum(Mat.new(1, col_size, self.vtype))
		when :row, :r
			_sum(Mat.new(row_size, 1, self.vtype))
		when :none, :n
			self.dup
		end
	end
	
	def mean(arg=:all)
		foo = sum(arg)
		if foo.kind_of?(Mat)
			if foo.row_size != self.row_size
				foo.e_div(self.row_size)
			else
				foo.e_div(self.col_size)
			end
		elsif arg.kind_of?(Mat)
			arg.e_div(self.length)
			arg[0, 0]
		else
			foo.quo(self.length)
		end
	end
	
	private def _mean(arg, ddof)
		foo = sum(arg)
		if foo.kind_of?(Mat)
			if foo.row_size != self.row_size
				foo.e_div(self.row_size-ddof)
			else
				foo.e_div(self.col_size-ddof)
			end
		elsif arg.kind_of?(Mat)
			arg.e_div(self.length-ddof)
			arg[0, 0]
		else
			foo.quo(self.length)
		end
	end
	def var(arg=:all, ddof: 0)
		foo = mean(arg)
		if foo.kind_of?(Mat)
			bar = foo.repmat(self.row_size.div(foo.col_size), self.col_size.div(foo.col_size))
			bar.sub!(self).e_mul!(bar)
		else
			unless arg.kind_of?(Mat)
				arg = Mat.new(1, 1, self.vtype) do
					foo
				end
			end
			bar = Mat.new(self.row_size, self.col_size, self.vtype)
			bar.fill(foo)
		end
		bar.sub!(self)
		bar.e_mul!(bar)
		bar.__send__(:_mean, arg, ddof)
		if foo.kind_of?(Mat)
			bar
		else
			bar[0, 0]
		end
	end
	alias :variance :var
	
	def std(arg=:all, ddof: 0)
		foo = var(arg, ddof: ddof)
		if foo.kind_of?(Mat)
			foo.sqrt!
		elsif arg.kind_of?(Mat)
			arg.sqrt!
			arg[0, 0]
		else
			Math.sqrt(arg)
		end
	end
	alias :stddev :std
	alias :standard_deviation :std
end
