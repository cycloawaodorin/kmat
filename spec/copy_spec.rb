RSpec.describe 'Copying elements:' do
	let(:a){ Mat.randn(2, 2) }
	let(:b){ Mat.randn(2, 2) }
	example '#dup' do
		ad = a.dup
		expect(ad).to eq a
		ad[0, 0] += 100
		expect(ad).not_to eq a
	end
	example '#copy_from' do
		a.copy_from(b)
		expect(a).to eq b
		expect(a.supermatrix).to be nil
	end
	example '#replace' do
		a.replace(b[])
		expect(a).to eq b
		expect(a.supermatrix).to be b
	end
	example '#fill' do
		a.fill(1.5)
		expect(a).to eq Mat[[1.5, 1.5], [1.5, 1.5]]
	end
	example 'marshal' do
		am = Marshal.load(Marshal.dump(a))
		expect(am).to eq a
		am[0, 0] += 100
		expect(am).not_to eq a
	end
end
