RSpec.describe 'Arithmetic operations:' do
	let(:a){ Mat.ones(2, 2) }
	let(:b){ Mat[[2, 0], [0, 1]] }
	example '#+' do
		expect(a+b).to eq Mat[[3, 1], [1, 2]]
		expect(b+2).to eq Mat[[4, 2], [2, 3]]
	end
	example '#-' do
		expect(a-b).to eq Mat[[-1, 1], [1, 0]]
		expect(b-1).to eq Mat[[1, -1], [-1, 0]]
	end
	example '#mprod' do
		expect(a.mprod(b)).to eq Mat[[2, 1], [2, 1]]
	end
	example '#e_mul' do
		expect(a.e_mul(b)).to eq Mat[[2, 0], [0, 1]]
	end
	example '#under' do
		expect(b.under(a)).to eq Mat[[0.5, 0.5], [1, 1]]
	end
	example '#over' do
		expect(a.over(b)).to eq Mat[[0.5, 1], [0.5, 1]]
	end
	example '#e_div' do
		expect{a.e_div(b)}.to raise_error(ZeroDivisionError)
		expect(b.e_div(a*2)).to eq Mat[[1, 0], [0, 0.5]]
	end
	example '#*' do
		expect{a*b}.to raise_error(ArgumentError)
		expect(b*2).to eq Mat[[4, 0], [0, 2]]
	end
	example '#/' do
		expect{a/b}.to raise_error(ArgumentError)
		expect(b/2).to eq Mat[[1, 0], [0, 0.5]]
	end
	context 'using Mat::MatrixProductOperator' do
		using Mat::MatrixProductOperator
		example '#*' do
			expect(a*b).to eq a.mprod(b)
		end
		example '#/' do
			expect(a/b).to eq a.over(b)
		end
	end
end
