import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from RNAPII_occupancy import STM, TFO

def Fit(model):

	#probe dimension of model output to see if it is TFO or STM
	test = model([1,1,1,1]).shape[0]

	params = []
	if test == 4: #STM
		#set up grid
		k4 = np.arange(0.0001,0.2002,0.01)
		k2 = np.arange(0.0001,0.2002,0.01)
		k2r = np.arange(0.0001,0.2002,0.01)
		k3r = np.arange(0.0001,0.0302,0.003)
		for ia in range(len(k4)):
		    for ib in range(len(k3r)):
		        for ic in range(len(k2)):
		            for id in range(len(k2r)):
		                params.append([k4[ia], k3r[ib], k2[ic], k2r[id]])

	if test == 3: #TFO
		#set up grid
		k4 = np.arange(0.0001,0.2002,0.001)
		k3r = np.arange(0.0001,0.0302,0.001)

		for ia in range(len(k4)):
		    for ib in range(len(k3r)):
		        params.append([k4[ia], k3r[ib]])

	params = np.array(params)


	#run over grid
	results = []
	for i in tqdm(range(len(params))):
		results.append(model(params[i,:]))
	results = np.array(results)


	#evaluate
	data = None
	if test == 3: #TFO
		data = np.array([0.4503438137077112,	0.048032317704999085,	0.20655627685074285]) #mean ChEC-seq2 data
	if test == 4: #STM
		data = np.array([0.34024874645180675,	0.46417911724809846,	0.066251271,	0.24924265793509923])#mean ChEC-seq2 data


	cosine_sim = []
	for i in range(len(params)):
		cosine_sim.append(np.dot(results[i,:], data)/np.linalg.norm(results[i,:])/np.linalg.norm(data))
	cosine_sim=np.array(cosine_sim)

	#threshold away poor fits
	thr = 0.995
	params = params[cosine_sim>0.995,:]


	#plot fit parameters
	if test == 4: #STM
		fig = plt.figure()
		ax = fig.add_subplot(projection='3d')
		im = ax.scatter(params[:,1], params[:,2], params[:,3], c=params[:,0], s=50)


		ax.set_xlabel(r"$K_{-3}$",size=12)
		ax.set_ylabel(r"$K_2$",size=12)
		ax.set_zlabel(r"$K_{-2}$",size=12)
		ax.set_xlim(0,0.03)
		ax.set_ylim(0,0.05)
		ax.set_zlim(0,0.05)

		ax.set_xticks([0,0.01, 0.02, 0.03])
		ax.set_yticks([0, 0.1,0.2])
		ax.set_zticks([0, 0.1,0.2])

		im.set_clim(0, 0.2)
		ax.view_init(elev=14., azim=-130)

		fig.tight_layout()
		plt.rc('font', size=14) 
		plt.show()

	if test == 3: #TFO
		fig = plt.figure()
		ax = fig.add_subplot()
		im = ax.scatter(params[:,0], params[:,1], s=50)

		ax.set_xlabel(r"$K_{4}$",size=12)
		ax.set_ylabel(r"$K_{-3}$",size=12)
		ax.set_ylim(0,0.03)
		ax.set_xlim(0,0.05)

		fig.tight_layout()
		plt.rc('font', size=14) 
		plt.show()

def main():
	print("Fitting STM model...")
	Fit(STM)

if __name__ == "__main__":
    main()









