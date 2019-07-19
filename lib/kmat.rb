require "kmat/version"

class Mat
	include Enumerable
end

class << Mat
	private def _define_dup_escaped_method(method)
		Mat.class_eval("def #{method}(*args, &block); self.dup.#{method}!(*args, &block); end")
	end
end

require "kmat/kmat"
require "kmat/arith"
require "kmat/accessor"
require "kmat/linalg"
require "kmat/logical"
require "kmat/misc"
require "kmat/random"
require "kmat/statistics"
