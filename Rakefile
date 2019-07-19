require "bundler/gem_tasks"
require "rake/extensiontask"

task :build => :compile

def istrue(key)
	val = ENV[key]
	if val.nil?
		false
	elsif val.kind_of?(String)
		if val == '0' || val =~ /false/i || val =~ /nil/i
			false
		else
			true
		end
	else
		raise "unknown environment value type #{val.class} for ENV['#{key}']"
	end
end

Rake::ExtensionTask.new("kmat") do |ext|
  ext.lib_dir = "lib/kmat"
  ext.config_options = %w|debug| if istrue('KMAT_DEBUG')
end

task :default => [:clobber, :compile, :spec]
