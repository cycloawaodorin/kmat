
class MethodDefinition
	def initialize(name_arg, comment)
		comment = %r|//\s+(.*)|.match(comment)[1]
		@als = /\Aalias (.+)/.match(comment)&.[](1)&.split(/[,\s]\s*/)
		@meth = /\Akmm_/.match(name_arg)
		return unless @meth
		m = %r|\A(kmm_([^\(_]+)_([^\(]+))\((.+)\)\Z|.match(name_arg)
		raise "unknown name_arg pattern `#{name_arg}' found" unless m
		@funcname, @type, @name, arg =  m[1], type_trans(m[2]), m[3], m[4]
		if m = %r|(.+)_p\Z|.match(@name)
			@name = m[1]+'?'
		elsif m = %r|(.+)_dest\Z|.match(@name)
			@name = m[1]
			@dup_esc = true
		elsif m = %r|(.+)_destl\Z|.match(@name)
			@name = m[1]+'!'
		end
		if m = %r|(.+)_m2\Z|.match(@name)
			@name = m[1]
			@argc = -2
		elsif m = %r|\Aint argc|.match(arg)
			@argc = -1
		else
			@argc = arg.split(/,/).size-1
		end
	end
	attr_reader :meth
	
	def type_trans(t)
		case t
		when 'obj'
			'rb_cObject'
		when 'mat'
			'km_cMat'
		when 'Mat'
			'km_sMat'
		when 'ary'
			'rb_cArray'
		when 'MATH'
			['rb_mMath', 'rb_sMath']
		when 'float'
			'rb_cFloat'
		else
			raise "unknown type `#{t}' found"
		end
	end
	
	def to_s(join="\n\t")
		ret = []
		if @type.kind_of?(Array)
			ret << %Q|rb_define_module_function(#{@type[0]}, "#{@name}", #{@funcname}, #{@argc});|
		elsif %r|\A\_|.match(@name)
			ret << %Q|rb_define_private_method(#{@type}, "#{@name}", #{@funcname}, #{@argc});|
		elsif @dup_esc
			ret << %Q|rb_define_method(#{@type}, "#{@name}!", #{@funcname}, #{@argc});|
			ret << %Q|rb_funcall(km_cMat, id__define_dup_escaped_method, 1, rb_str_new_cstr("#{@name}"));|
		else
			ret << %Q|rb_define_method(#{@type}, "#{@name}", #{@funcname}, #{@argc});|
		end
		@als&.each do |a|
			if @type.kind_of?(Array)
				ret << %Q|rb_define_alias(#{@type[0]}, "#{a}", "#{@name}");|
				ret << %Q|rb_define_alias(#{@type[1]}, "#{a}", "#{@name}");|
			else
				ret << %Q|rb_define_alias(#{@type}, "#{a}", "#{@name}");|
				ret << %Q|rb_define_alias(#{@type}, "#{a}!", "#{@name}!");| if @dup_esc
			end
		end
		ret.join(join)
	end
end

defs, decs = [], []

# for all .c files
Dir.glob("#{__dir__}/**/*.c").each do |file|
	# for lines that comment line (can be omitted), return type and modifier line, function name and argument line
	File.read(file).gsub(%r|(//[^\n]+\n)?[^#/{}\n\s][^\n;]+\n[^#/{}\n\s][^\n;]+\n{\n|) do |m|
		ary = m.split(/\n/)
		# functions which name start with kmm_ are methods
		if ary.size == 4
			foo = MethodDefinition.new(ary[2], ary[0])
		else
			foo = MethodDefinition.new(ary[1], '// ')
		end
		defs << foo.to_s if foo.meth
		# ignore static functions and static variables
		unless /\Astatic/.match(ary[-3])
			decs << "#{ary[-3]} #{ary[-2]};"
		end
		''
	end
end


File.open('./auto_collected.h', 'w') do |f|
	decs.each do |dec|
		f.puts dec
	end
end

File.open('./method_definitions.c', 'w') do |f|
	f.puts "static void\nkm_define_methods(void)\n{"
	defs.each do |def_|
		f.puts "\t#{def_}"
	end
	f.puts '}'
end

File.open('./global_variables.h', 'w') do |f|
	flg = false
	File.foreach("#{__dir__}/main.c") do |line|
		flg = true if line == "// km_global_variables_begin\n"
		f.puts "extern #{line}" if flg && !%r|\A// km_global_variables|.match(line) && !%r|\A\s*\n\Z|.match(line)
		break if line == "// km_global_variables_end\n"
	end
end
