#include <fstream>

void Getevtnb(const char* dataType = "", bool pnfs = true){
std::string chain_name = "bTag_AntiKt4HIJets";
   TTree *myChain = (TTree*)f->Get(chain_name.c_str());

   if (pnfs){
	   ifstream
   }
	DIR *dir1;
	dirent *pdir;
	dir1 = opendir(input.c_str());
	while ((pdir = readdir(dir1)))
	{
		std::string foldName = pdir->d_name;
		cout << pdir->d_name << endl;
		if (foldName.find(inputFolder) == std::string::npos)
			continue;
		if (!(foldName.find(".root") == std::string::npos))
			continue;
		//if(foldName.find("JZ5")==std::string::npos) continue;
		cout << "Success:" << pdir->d_name << endl;
		DIR *dir2;
		dirent *pdir2;
		dir2 = opendir((input + "/" + foldName).c_str());
		while ((pdir2 = readdir(dir2)))
		{
			std::string fName = pdir2->d_name;
			if (fName.find(".root") == std::string::npos)
				continue;

			//if (fName.find("JZ5") == std::string::npos)
			//	continue;
			myChain->Add((input + "/" + foldName + "/" + fName).c_str());
			//goto here;
		}
	}
 
