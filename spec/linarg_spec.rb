RSpec.describe 'Linear algebraic operations:' do
	let(:m){ 3 }
	let(:n){ 3 }
	let(:a){ Mat.randn(m, n) }
	let(:b){ Mat.rand(n, 1) }
	let(:eps){ 1024*Float::EPSILON }
	let(:eye){ Mat.eye(m, n) }
	let(:zeros){ Mat.zeros(m, n) }
	example '#ls' do
		x = a.ls(b)
		expect(a.mprod(x).near?(b, eps)).to be true
	end
	example '#evd' do
		v, d = a.evd
		expect(a.to_cmat.mprod(v).near?(v.mprod(d), eps)).to be true
		d2 = a.eigen_values
		expect(d.diag.near?(d2, eps))
		d.diag[] = 0
		expect(d).to eq zeros.to_cmat
	end
	example '#sym_evd' do
		a_sym = a.symmetrize
		v, d = a_sym.evd
		expect(v.mprod(d).mprod(v.t).near?(a_sym, eps)).to be true
		expect(v.mprod(v.t).near?(eye, eps)).to be true
		d.diag[] = 0
		expect(d).to eq zeros
	end
	example '#svd' do
		u, s, v = a.svd
		expect(a.near?(u.mprod(s).mprod(v.t), eps)).to be true
		expect(u.mprod(u.t).near?(eye, eps)).to be true
		expect(v.mprod(v.t).near?(eye, eps)).to be true
		s2 = a.singular_values
		expect(s.diag.near?(s2, eps)).to be true
		s.diag[] = 0
		expect(s).to eq zeros
	end
	example '#lu' do
		l, u = a.lu
		expect(a.near?(l.mprod(u), eps)).to be true
		(-1).downto(-m+1) do |k|
			expect(u.diag(k)).to eq Mat.zeros(m+k, 1)
		end
		expect(l.count{|e| e==0.0}).to be >= 3
	end
	example '#lup' do
		l, u, pm = a.lup
		expect(pm.mprod(a).near?(l.mprod(u), eps)).to be true
		1.upto(n-1) do |k|
			expect(l.diag(k)).to eq Mat.zeros(n-k, 1)
		end
		(-1).downto(-m+1) do |k|
			expect(u.diag(k)).to eq Mat.zeros(m+k, 1)
		end
		expect(pm.count{|e| e!=0.0 && e!=1.0}).to be 0
		expect(pm.sum(:row)).to eq Mat.ones(m, 1)
		expect(pm.sum(:col)).to eq Mat.ones(1, n)
	end
	example '#det' do
		expect(a.det.abs).to be_within(eps).of(a.singular_values.prod)
	end
	example '#qr' do
		q, r = a.qr
		expect(a.near?(q.mprod(r), eps)).to be true
		expect(q.mprod(q.t).near?(eye, eps)).to be true
		(-1).downto(-m+1) do |k|
			expect(r.diag(k)).to eq Mat.zeros(m+k, 1)
		end
	end
	example '#expm' do
		a_sym = a.symmetrize
		v, d = a_sym.sym_evd
		expect(a_sym.expm.near?(v.mprod(Mat.diag(d.diag.exp)).mprod(v.t), eps)).to be true
	end
	example '#chol' do
		a_sym = a.mprod(a.t).symmetrize
		u = a_sym.chol
		expect(a_sym.near?(u.t.mprod(u), eps)).to be true
		(-1).downto(-m+1) do |k|
			expect(u.diag(k)).to eq Mat.zeros(m+k, 1)
		end
	end
end
