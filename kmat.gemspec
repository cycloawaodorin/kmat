
lib = File.expand_path("../lib", __FILE__)
$LOAD_PATH.unshift(lib) unless $LOAD_PATH.include?(lib)
require "kmat/version"

Gem::Specification.new do |spec|
  spec.name          = "kmat"
  spec.version       = Mat::VERSION
  spec.authors       = ["KAZOON"]
  spec.email         = ["cycloawaodorin+gem@gmail.com"]
  spec.license = "GPL-3.0"

  spec.summary       = %q{Kmat is a Ruby gem for matrix operations. Kmat uses BLAS/LAPACK as back-end.}
  #spec.description   = %q{TODO: Write a longer description or delete this line.}
  spec.homepage      = "https://github.com/cycloawaodorin/"

  # Prevent pushing this gem to RubyGems.org. To allow pushes either set the 'allowed_push_host'
  # to allow pushing to a single host or delete this section to allow pushing to any host.
  if spec.respond_to?(:metadata)
    #spec.metadata["allowed_push_host"] = "http://localhost"

    spec.metadata["homepage_uri"] = spec.homepage
    spec.metadata["source_code_uri"] = "https://github.com/cycloawaodorin/kmat"
    spec.metadata["changelog_uri"] = "https://github.com/cycloawaodorin/kmat/blob/master/CHANGELOG.md"
  else
    raise "RubyGems 2.0 or newer is required to protect against " \
      "public gem pushes."
  end

  # Specify which files should be added to the gem when it is released.
  # The `git ls-files -z` loads the files in the RubyGem that have been added into git.
  spec.files         = Dir.chdir(File.expand_path('..', __FILE__)) do
    `git ls-files -z`.split("\x0").reject { |f| f.match(%r{^(test|spec|features)/}) }
  end
  spec.bindir        = "exe"
  spec.executables   = spec.files.grep(%r{^exe/}) { |f| File.basename(f) }
  spec.require_paths = ["lib"]
  spec.extensions    = ["ext/kmat/extconf.rb"]

  spec.add_development_dependency "bundler", ">= 2.4.8"
  spec.add_development_dependency "rake", ">= 13.0.6"
  spec.add_development_dependency "rake-compiler"
  spec.add_development_dependency "pry", ">= 0.14.2"
  
  spec.required_ruby_version = '>= 3.2.0'
end
