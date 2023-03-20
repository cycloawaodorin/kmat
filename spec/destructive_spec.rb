RSpec.describe 'Destructive operations:' do
	let(:a){ Mat.ones(2, 2) }
	let(:b){ Mat[[2, 0], [0, 1]] }
	let(:c){ Mat.new(2, 2) }
	let(:a_org){ Mat.ones(2, 2) }
	let(:b_org){ Mat[[2, 0], [0, 1]] }
	example '#add!' do
		_ = a+b
		expect([a, b]).to eq [a_org, b_org]
		a.add!(b)
		expect([a, b]).to eq [Mat[[3, 1], [1, 2]], b_org]
	end
	example '#sub!' do
		_ = a-b
		expect([a, b]).to eq [a_org, b_org]
		a.sub!(b)
		expect([a, b]).to eq [Mat[[-1, 1], [1, 0]], b_org]
	end
	example '#e_mul!' do
		a.e_mul(b)
		expect([a, b]).to eq [a_org, b_org]
		a.e_mul!(b)
		expect([a, b]).to eq [Mat[[2, 0], [0, 1]], b_org]
	end
	example '#e_div!' do
		b.e_div(a*2)
		expect([a, b]).to eq [a_org, b_org]
		b.e_div!(a*2)
		expect([a, b]).to eq [a_org, Mat[[1, 0], [0, 0.5]]]
	end
	example '#mprod!' do
		a.mprod(b)
		expect([a, b]).to eq [a_org, b_org]
		c.mprod!(a, b)
		expect([a, b, c]).to eq [a_org, b_org, Mat[[2, 1], [2, 1]]]
	end
	example '#under!' do
		b.under(a)
		expect([a, b]).to eq [a_org, b_org]
		c.under!(b, a)
		expect([a, b, c]).to eq [a_org, b_org, Mat[[0.5, 0.5], [1, 1]]]
	end
	example '#over!' do
		a.over(b)
		expect([a, b]).to eq [a_org, b_org]
		c.over!(a, b)
		expect([a, b, c]).to eq [a_org, b_org, Mat[[0.5, 1], [0.5, 1]]]
	end
	example '#s_mul!' do
		_ = b*2
		expect(b).to eq b_org
		b.s_mul!(2)
		expect(b).to eq Mat[[4, 0], [0, 2]]
	end
	example '#s_div!' do
		_ = b/2
		expect(b).to eq b_org
		b.s_div!(2)
		expect(b).to eq Mat[[1, 0], [0, 0.5]]
	end
end
