RSpec.describe 'Sort, stacking and logical operations:' do
	let(:a){ Mat[[3, 2, 1], [5, -3, 7]] }
	let(:b){ Mat[[1], [5]] }
	example 'sort and flip' do
		expect(a.sort(:row)).to eq Mat[[1, 2, 3], [-3, 5, 7]]
		expect(a.rsort(:col)).to eq Mat[[5, 2, 7], [3, -3, 1]]
		expect(a.flip(:both)).to eq Mat[[7, -3, 5], [1, 2, 3]]
		expect(a.t).to eq Mat[[3, 5], [2, -3], [1, 7]]
	end
	example 'stacking' do
		expect(Mat.blocks([[a, b], [b, a]])).to eq Mat[[3, 2, 1, 1], [5, -3, 7, 5], [1, 3, 2, 1], [5, 5, -3, 7]]
		expect(Mat.vstack(a, a)).to eq Mat[[3, 2, 1], [5, -3, 7], [3, 2, 1], [5, -3, 7]]
		expect(Mat.hstack(a, b)).to eq Mat[[3, 2, 1, 1], [5, -3, 7, 5]]
	end
	example 'logical operations' do
		expect(a.gt(b).to_ary).to eq [[true, true, false], [false, false, true]]
		expect(a.eq(b).to_ary).to eq [[false, false, true], [true, false, false]]
		expect(a.ge(b).to_ary).to eq [[true, true, true], [true, false, true]]
		expect(a.lt(b).to_ary).to eq [[false, false, false], [false, true, false]]
		expect(a.le(b).to_ary).to eq [[false, false, true], [true, true, false]]
		expect(a.max).to eq 7.0
		expect(a.maximum(b)).to eq Mat[[3, 2, 1], [5, 5, 7]]
		expect{Mat.maximum(a, b, Mat.randn(2, 3, random: Random.new(0)))}.to raise_error Mat::MismatchedDimensionError
		expect(Mat.maximum(a, b.repmat(1, 3), Mat.randn(2, 3, random: Random.new(3))).near?(Mat[[3, 2, 1.6808342571962616], [5, 5, 7]])).to be true
	end
end
