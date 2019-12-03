class << Mat
	def hilb(n, type=:float)
		self.new(n, n, type) do |i, j|
			1.quo(i+j+1)
		end
	end
	def load(file)
		case file
		when String
			ret = nil
			File.open(file) do |f|
				ret = Marshal.load(f)
			end
			ret
		when IO
			Marshal.load(file)
		else
			raise ArgumentError, 'the argument must be a filepath or an IO object'
		end
	end
end

class Mat
	def save(file)
		case file
		when String
			File.open(file, 'w') do |f|
				Marshal.dump(self, f)
			end
		when IO
			Marshal.dump(self, file)
		else
			raise ArgumentError, 'the argument must be a filepath or an IO object'
		end
		nil
	end
	
	def geo_normalize!(arg=:all)
		self.log!
		foo = self.mean(arg)
		if foo.kind_of?(Mat)
			self.sub!(self.broadcast(foo))
		else
			self.s_sub!()
		end
		self.exp!
	end
	def geo_normalize(arg=:all)
		self.dup.geo_normalize!(arg)
	end
	def normalize!
		if self.vector?
			self.s_div!(self.normf)
		elsif self.square?
			ev = self.eigen_values
			if ev[0] < 0
				self.diag.s_add!(ev[0]*(-2))
				ev.s_add!(ev[0]*(-2))
			end
			self.s_div!(Math.exp(ev.log!.mean))
		else
			raise MismatchedDimensionError, "normalize is available only for vectors or square matricies"
		end
	end
	def normalize
		self.dup.normalize!
	end
end

class Array
	def argsort
		ret = Array.new(self.size, &:itself)
		if block_given?
			ret.sort! do |a, b|
				yield(self[a], self[b])
			end
		else
			ret.sort_by! do |elm|
				self[elm]
			end
		end
	end
	def argsort!(ary)
		if block_given?
			sort! do |a, b|
				yield(ary[a], ary[b])
			end
		else
			sort_by! do |elm|
				self[elm]
			end
		end
	end
	def argsort_by
		if block_given?
			Array.new(self.size, &:itself).sort_by! do |elm|
				yield(self[elm])
			end
		else
			self.to_enum(:argsort_by)
		end
	end
end
module Enumerable
	def argsort(&block)
		self.to_a.argsort(&block)
	end
	def argsort_by(&block)
		self.to_a.argsort_by(&block)
	end
end
class Numeric
	def limit(min, max)
		if ( self < min ) then
			min
		elsif ( max < self ) then
			max
		else
			self
		end
	end
end
