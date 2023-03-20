RSpec.describe 'Construction of matricies:' do
	example '.[]' do
		mat = Mat[[1, 2], [3, 4]]
		expect(mat).to be_an_instance_of Mat
		expect(mat.to_ary).to eq [[1.0, 2.0], [3.0, 4.0]]
	end
	example '.zeros' do
		expect(Mat.zeros(3, 2)).to eq Mat[[0, 0], [0, 0], [0, 0]]
	end
	example '.ones' do
		expect(Mat.ones(2, 3)).to eq Mat[[1, 1, 1], [1, 1, 1]]
	end
	example '.randn' do
		expect(Mat.randn(2, 2, random: Random.new(0)).near?(Mat[[-0.27375423029655194, -1.3051634279785924], [-1.2315874462339091, -0.37814642385629815]])).to be true
	end
	example '.range' do
		expect(Mat.range(3)).to eq Mat[[0], [1], [2]]
	end
	example '.new' do
		expect(Mat.new(2, 3){|i, j| i+j}).to eq Mat[[0, 1, 2], [1, 2, 3]]
	end
end
