require "mkmf"
require 'fileutils'

if find_library("mkl_rt", nil, "/opt/intel/mkl/lib/intel64")
	puts "It is set to use Intel MKL for BLAS/LAPACK functions."
	mkl = true
elsif have_library("blas") && have_library("lapack")
	puts "It is set to use BLAS/LAPACK."
	mkl = false
else
	puts "Intel MKL or BLAS/LAPACK is neaded to build this library."
	exit
end

$CFLAGS = "$(cflags) -std=c11"
$CFLAGS += " -m64" if mkl

#$warnflags = "-Wall -Wextra -Wdeprecated-declarations -Wimplicit-function-declaration -Wimplicit-int -Wpointer-arith -Wwrite-strings -Wmissing-noreturn -Wno-unused-parameter -Wsuggest-attribute=format -Wsuggest-attribute=noreturn -Wunused-variable -Wno-maybe-uninitialized -Winit-self -Wshadow"
$warnflags = "-Wall -Wextra -Wdeprecated-declarations -Wimplicit-function-declaration -Wimplicit-int -Wpointer-arith -Wwrite-strings -Wmissing-noreturn -Wno-unused-parameter -Wsuggest-attribute=format -Wsuggest-attribute=noreturn -Wunused-variable -Winit-self -Wshadow -Wlogical-op -Wconversion"

$DLDFLAGS += " -Wl,-Bsymbolic -fPIC"
$DLDFLAGS += " -Wl,--no-as-needed" if mkl


# set .c files in subdirectories as source
$objs = Dir.glob("#{__dir__}/**/*.c").map do |file|
	file[-1] = 'o'
	obj = %r|#{__dir__}/(.+)|.match(file)[1]
	if m = %r|(.+)/[^/]+|.match(obj)
		FileUtils.mkdir_p(m[1]) unless FileTest.exist?(m[1])
	end
	obj
end

srcs = $objs.map do |file|
	file = file.dup
	file[-1] = 'c'
	"$(srcdir)/#{file}"
end

create_makefile("kmat/kmat")

# change variables after Makefile created
File.open('./__Makefile__temp__', 'w') do |f|
	File.foreach('Makefile') do |line|
		if ARGV.include?('debug')
			line.sub!(/^optflags.+/, 'optflags = -O0')
		else
			line.sub!(/^debugflags.+/, 'debugflags =')
		end
		f.puts line
	end
end
File.unlink('./Makefile')
FileUtils.move('./__Makefile__temp__', './Makefile')

# add dependencies
File.open('./Makefile', 'a') do |f|
	f.puts "true_srcs = #{srcs.join(' ')}"

	# invoke auto_collect.rb
	f.puts "auto_collected.h: auto_collect.rb $(true_srcs)"
	f.puts "\truby $(srcdir)/auto_collect.rb"
	
	# invoke id_sym.rb
	f.puts "id_sym.c id_sym.h: id_sym.rb id.txt sym.txt"
	f.puts "\truby $(srcdir)/id_sym.rb"
	
	# invoke elementwise_functions.rb
	f.puts "elementwise_function.h: elementwise_function.rb"
	f.puts "\truby $(srcdir)/elementwise_function.rb"

	f.puts '$(OBJS): auto_collected.h id_sym.h elementwise_function.h'
	f.puts 'main.o: id_sym.c'
end
