RSpec.describe 'Value types:' do
	example 'float' do
		a = Mat.new(2, 2, :float){|i, j| i+j}
		expect([a.vtype, a.to_ary]).to eq [:float, [[0.0, 1.0], [1.0, 2.0]]]
		expect(a[0, 0]).to be_an_instance_of Float
	end
	example 'complex' do
		a = Mat.new(2, 2, :complex){|i, j| Complex.rect(i, j)}
		expect([a.vtype, a.to_ary]).to eq [:complex, [[0.0, 1.0i], [1.0, 1.0+1.0i]]]
		expect(a[0, 0]).to be_an_instance_of Complex
	end
	example 'int' do
		a = Mat.new(2, 2, :int){|i, j| i-j}
		expect([a.vtype, a.to_ary]).to eq [:int, [[0, -1], [1, 0]]]
		expect(a[0, 0]).to be_an_instance_of Integer
	end
	example 'bool' do
		a = Mat.new(2, 2, :bool){|i, j| i==j}
		expect([a.vtype, a.to_ary]).to eq [:bool, [[true, false], [false, true]]]
		expect(a[0, 0]).to be_an_instance_of TrueClass
		expect(a[0, 1]).to be_an_instance_of FalseClass
	end
	example 'object' do
		a = Mat.new(2, 2, :object){|i, j| Rational(i, j+1)}
		expect([a.vtype, a.to_ary]).to eq [:object, [[0, 0], [1, 1.quo(2)]]]
		expect(a[0, 0]).to be_an_instance_of Rational
	end
end
