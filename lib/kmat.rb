require_relative "./kmat/version"

class Mat
	include Enumerable
end

class << Mat
	private def _define_dup_escaped_method(method)
		Mat.class_eval("def #{method}(*args, &block); self.dup.#{method}!(*args, &block); end")
	end
end

require_relative "./kmat/kmat"
require_relative "./kmat/arith"
require_relative "./kmat/accessor"
require_relative "./kmat/linalg"
require_relative "./kmat/logical"
require_relative "./kmat/misc"
require_relative "./kmat/random"
require_relative "./kmat/statistics"
