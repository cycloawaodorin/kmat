RSpec.describe 'Numpy-like broadcasting:' do
	let(:a){ Mat.eye(2, 2) }
	let(:b){ Mat.range(2) }
	example 'non-destructive' do
		expect(a+b).to eq Mat[[1, 0], [1, 2]]
	end
	example 'destructive' do
		expect{a.add!(b)}.to raise_error Mat::MismatchedDimensionError
		a.add!(a.broadcast(b))
		expect(a).to eq Mat[[1, 0], [1, 2]]
	end
	example 'scalar product' do
		expect(0.5*a).to eq Mat[[0.5, 0], [0, 0.5]]
	end
end
