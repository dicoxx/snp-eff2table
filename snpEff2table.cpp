
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <math.h>
#include <cassert>
#include <unistd.h>
#include <getopt.h>
using namespace std;


void ShowUsage()
{
	cout << "\nUsage: snpEff2table <OPTION>\n" << endl;
	cout << "Option:\n" << endl;
	cout << "\t-i \tinput (snpEff output)" << endl;
	cout << "\t-o \toutput "<< endl;
	cout << " "<< endl;
	return ;
}


using namespace std;


string simpleSub(string simple)
{
	string str = ":";
	int lenBegin_6 = simple.find_last_of(str);
	string begin_6_1 = simple.substr(0, lenBegin_6);
	int lenBegin_6_1 = begin_6_1.find_last_of(str);
	string begin_6_2 = begin_6_1.substr(0, lenBegin_6_1); 
	return begin_6_2;

}
vector<string>sub_Sev;  //存放每一行第七列分割后的数据

void split(string simple,string delimiter,vector<string>&sub_Sev)
{
	int posBegin = simple.find(delimiter);
	if (posBegin != simple.npos) {
		string tmp = simple.substr(0, posBegin);
		string last = simple.substr(posBegin+1);
		sub_Sev.push_back(tmp);
		split(last, delimiter,sub_Sev);

	}
	else
	{
		return ;
	}
}
void snpEff2table(string srcFile,string dstFile)
{
	ifstream inFile;
	inFile.open(srcFile, ios::in | ios::out);
	string str;
	vector<vector<string>>in;
	while (getline(inFile, str))
	{

		istringstream input(str);
		vector<string> tmp;
		string tmp_s;
		while (input >> tmp_s)
		{
			tmp.push_back(tmp_s);
		}
		in.push_back(tmp);
	}
	/*for(int i=0;i<97;i++)
	{
		for(int j=0;j<in[i].size();j++)
		{
			cout<<in[i][j]<<"\t";
		}
		cout<<endl;
	}*/
	vector < vector<string>>out;
	int tag;
	for (int i = 0; i< in.size(); i++)
	{
		for (int j = 0; j < in[i].size(); j++)
		{
			int a=in[i][j].find("POS");
			if (a != -1)
				tag = i;
		}

	}
	vector<string> name;
	name.push_back("Vars");
	name.push_back("CHROM");
	name.push_back("POS");
	name.push_back("REF");
	name.push_back("ALT");
	name.push_back("Effect");
	name.push_back("Position");
	name.push_back("Transcript");
	for (int j = 9; j < in[tag].size(); j++)
	{
		name.push_back(in[tag][j] + ".GT");
		name.push_back(in[tag][j] + ".AD");
		name.push_back(in[tag][j] + ".DP");
	}

	out.push_back(name);

	vector<vector<string>>Seventh;  //存放所有数据第七列分割后的数据
	/*for (int i = 0; i < in.size(); i++)
	{
		vector<string>p;

		if(in[i][0].substr(0,3) != "Chr")
		{
			continue;
		}
		else
		{
			split(in[i][7], "|", sub_Sev);
			Seventh.push_back(sub_Sev);
		}
		sub_Sev.clear();
	}*/
	int tag_i=0;
	for (int i = 0; i < in.size(); i++)
	{
		if(in[i][0].substr(0,6) == "#CHROM")
		{
			tag_i=i;
			break;
		}
	}
	for (int i = tag_i+1; i < in.size(); i++)
	{
		vector<string>p;

		//if(in[i][0].substr(0,3) == "Chr" || in[i][0].substr(0,1) == "A")
		{
			split(in[i][7], "|", sub_Sev);
			Seventh.push_back(sub_Sev);
		}
		sub_Sev.clear();
	}
	/*for(int i=0;i<Seventh.size();i++)
	{
		for(int j=0;j<Seventh[i].size();j++)
		{
			cout<<Seventh[i][j]<<"\t";
		}
		cout<<endl;
	}*/

	for (int i = 0; i < in.size(); i++)
	{

		vector<string> vv;
		for (int j = 0; j < in[i].size(); j++)
		{
			/*if (in[i][0].substr(0, 3) == "Chr"||in[i][0].substr(0,1) == "A")
			{
				break;
			}
			else*/
			if(i>tag_i)
			{

				//第一列
				vv.push_back(in[i][0] + "_" + in[i][1] + "_" + in[i][3] + "_" + in[i][4]);
				//CHROM
				vv.push_back(in[i][0]);
				//POS
				vv.push_back(in[i][1]);
				//REF
				vv.push_back(in[i][3]);
				//ALT
				vv.push_back(in[i][4]);
				//Effect
				vv.push_back(Seventh[i-tag-1][1]);
				//Position
				vv.push_back(Seventh[i-tag-1][9]);
				//transcript
				vv.push_back(Seventh[i-tag-1][6]);

				//样本数据列
				for (j = 9; j < in[i].size();j++)
				{
					string simple;
					simple=simpleSub(in[i][j]);
					//simple.DP
					string str3 = simple.substr(simple.find_last_of(":")+1);
					//simple.GT
					string tmp1 = simple.substr(0, simple.find_last_of(":"));
					string str1 = tmp1.substr(0, tmp1.find(":"));
					str1 = str1.replace(str1.find("/"), 1, "|");
					//simple.AD
					string str2 = simple.substr(simple.find(":")+1, simple.find_last_of(":") - simple.find(":")-1);
					str2 = str2.replace(str2.find_first_not_of("0123456789"), 1, "|");
					vv.push_back(str1);
					vv.push_back(str2);
					vv.push_back(str3);

				}
			//cout<<"第"<<i<<"行尺寸为：	"<<vv.size()<<endl;
				//cout<<vv[i]<<endl;
				out.push_back(vv);
				vector<string>().swap(vv);
				break;


			}
		}
	}


	std::ofstream outFile(dstFile, ios::out);
	/*for (int i = 0; i < out.size(); i++)
	{

		for (int j = 0; j < out[i].size(); j++)
		{
			cout << out[i][j] << "\t";
		}
		cout << endl;
	}*/


	for (int i = 0; i<out.size(); i++)
	{
		for (int j = 0; j<out[i].size(); j++)
		{
			if(j==out[i].size()-1)
			{
				outFile << out[i][j];
			}else
			{
				outFile << out[i][j] << "\t";
			}
		}
		outFile << endl;
	}

	inFile.close();
        outFile.close();

	return ;


}

int main(int argc, char* argv[])
{
	if(argc < 4)
	{
		ShowUsage();
		cout<<"\tExample: snpEff2table -i input.txt -o output.txt"<<endl;
		cout<<""<<endl;
		return -1;
	}
	int lopt;

	static struct option long_options[]={
		{"input",required_argument,NULL,'i'},
		{"output",required_argument,NULL,'o'},
		{0,0,0,0}
	};
	static char* const short_options=(char *)"i:o:";
	int option_index=0;

	string strInput;  // 存储输入文件名
	string strOutput; // 临时存储输入文件前缀，筛选方式以及筛选临界阈值
	int in=0;
	while((in=getopt_long(argc,argv,short_options,long_options,&option_index)) != -1)
	{

		ostringstream oss;
		switch(in){
			case 'i':
			       strInput=optarg;
			       break;
			case 'o':
			       strOutput=optarg;
			       break;
			case '?':
			       cout<<"unknown command"<<endl;
			       return -1;
			       //break;
		}
	}
	snpEff2table(strInput,strOutput);
	return 0;
};
