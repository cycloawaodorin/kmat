RSpec.describe 'Element or submatrix accessing:' do
	let(:a){ Mat[[1, 2], [3, 4]] }
	let(:vidx){ Mat.new(1, 2, :object){ |i, j| [j, i] } }
	example 'element accessing' do
		expect(a[1, 0]).to eq 3.0
		expect(a[0, -1]).to eq 2.0
		expect(a[nil, 1]).to eq Mat[[2], [4]]
		expect(a[[1, 0], 0..1]).to eq Mat[[3, 4], [1, 2]]
		expect(a[a.gt(2)]).to eq Mat[[3], [4]]
		expect(a[vidx]).to eq Mat[[1, 3]]
	end
	example 'submatrix making' do
		b = a[0, nil]
		expect(b).to eq Mat[[1, 2]]
		a[0, 0] = 5
		expect(b).to eq Mat[[5, 2]]
		b[0, 1] = -1
		expect(a).to eq Mat[[5, -1], [3, 4]]
	end
	example 'freeze' do
		b = a[0, nil]
		b.freeze
		expect{a[0, 0] = -3}.not_to raise_error
		a.deep_freeze
		expect{b[0, 0] = 7}.to raise_error FrozenError
		c = a[0, nil]
		expect{c[0, 0] = 9}.to raise_error FrozenError
	end
	example 'other types of accessing' do
		a[nil, 0] = Mat.range(2)
		expect(a).to eq Mat[[0, 2], [1, 4]]
		expect(a.map(&:-@)).to eq Mat[[0, -2], [-1, -4]]
		a.map_with_index!{|e, i, j| e+j}
		expect(a.diag).to eq Mat[[0], [5]]
		expect(a[]).to eq Mat[[0, 3], [1, 5]]
		expect(Mat.range(3)[2]).to eq 2.0
	end
end
