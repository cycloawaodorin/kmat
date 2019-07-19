id_list, sym_list = Hash.new, Hash.new
[['id', id_list], ['sym', sym_list]].each do |is, list|
	File.foreach("#{__dir__}/#{is}.txt").with_index do |line, i|
		ary = line.chomp.split(' ')
		case ary.size
		when 0
			nil
		when 1
			list.store(ary[0], ary[0])
		when 2
			list.store(ary[0], ary[1])
		else
			raise "Line #{i+1} in #{is}.txt cannot be recognized"
		end
	end
end

File.open('id_sym.c', 'w') do |f|
	[['id', id_list, 'ID', ''], ['sym', sym_list, 'VALUE', 'rb_id2sym']].each do |is, list, iv, func|
		list.each_key do |key|
			f.puts "#{iv} #{is}_#{key};"
		end
		f.puts ''
		f.puts "void\nkm_init_#{is}(void)\n{"
		list.each do |key, value|
			f.puts "\t#{is}_#{key} = #{func}(rb_intern(#{value.inspect}));"
		end
		f.puts '}'
		f.puts ''
	end
end

File.open('id_sym.h', 'w') do |f|
	[['id', id_list, 'ID'], ['sym', sym_list, 'VALUE']].each do |is, list, iv|
		list.each_key do |key|
			f.puts "extern #{iv} #{is}_#{key};"
		end
		f.puts ''
	end
end
