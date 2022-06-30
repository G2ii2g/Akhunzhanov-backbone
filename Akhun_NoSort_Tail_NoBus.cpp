#define _CRT_SECURE_NO_WARNINGS
#include <iostream>
#include <string>
#include <math.h>
#include <vector>
#include <fstream>
#include <sstream>
#include <ctime>
#include <algorithm>
#include <profileapi.h>

using namespace std;
const int black = -1;
const int green = -2;
const int yellow = -3;
// ------------------------
// ����� �������� � ��������� ������ - ����� ��������
char* inputFileName = (char*)"graph.txt";
char* resultFileName = (char*)"result.txt";
// ------------------------
#ifdef DEBUG
ofstream fout ("screen.txt");
#endif // DEBUG
struct Vertex {
    int row;
    int col;
    int color;
    vector<int> edges;
    Vertex(int i, int j, int c) {
        row = i;
        col = j;
        color = c;
    }
};
struct Edge {
    int v1;
    int v2;
    int color;

    Edge(int i1, int i2, int c) {
        v1 = i1;
        v2 = i2;
        color = c;
    }
};
//-----------------------------------------------------------
vector<Vertex> vertexList;
vector<Edge> edgeList;

static inline int64_t GetTicks()
{
    LARGE_INTEGER ticks;
    QueryPerformanceCounter(&ticks);

    return ticks.QuadPart;
}

static inline int64_t GetTicksPerSec()
{
    LARGE_INTEGER Frequency;

    QueryPerformanceFrequency(&Frequency);
    return Frequency.QuadPart;
}

int NextEdge(vector<int> &edges, int iRel) {// edges - ������ ����� �������, iRel- ������ ����� � ������ ����� ������������� � �������
    //������� ���������� ������ � edgeList ���������� ����� �� ������� �������, ���������� ��������� ����� �� ����� � ������������,
    //��� ����� � ������ �����, ������������� � �������, ������������� ������ ������� �������
    if (iRel > 0) return  edges[iRel - 1];
    else return edges[edges.size() - 1];
}

int otherVertex(int iV, int iE) {// ������� � �������� iV, ����������� ����� � �������� iE.
    // ������� ���������� ������ ������� ����������� ����� � �������� iE
    Edge& start = edgeList[iE];
    int iv1 = start.v1;//  ������ ������� ����� ����� iE �� ����������� � iV
    if (iV == iv1) iv1 = start.v2;
    return iv1;
}
#ifdef DEBUG
void printvertex(int i) {
    Vertex& V = vertexList[i];
    fout << i << " " << V.row << " " << V.col << " " << V.color << " " << V.edges.size() << ":";
    for (size_t j = 0; j < V.edges.size(); j++) {
        fout << V.edges[j] << ",";
    }
    fout << endl;
}
#endif // DEBUG
int IndexOfEdge(vector<int> &edges , int iE){//������� ���������� ������ ����� iE � ������ ����� ������� iV
    // ���� ������ ����� ���, �� ������� ���������� -1
    int index = -1;
    for (size_t i=0;i<edges.size();i++)
        if (edges[i] == iE) {
            index = i;
            break;
        }
    return index;
}

int NextYellowEdge(vector<int> &edges, int index){// ������ ������ ������ ����� ����� ����� � �������� index � ������ edges
    //����� ������� �� ������� �������, ���� ������ ����� ����� ������� index ���, �� ������������ index,
    //���� index ����������� � ������ �����, �� ��������� �����������
    if(index<0)
        exit (1);//������������ ��������: ����� ����������� � ������ �����
    int length = edges.size();
    for (int i = length -1; i > 0;i--) {// ������ ������ ������ ����� �� ����� � �������� index �� ������� �������,
    //���� ������ ����� ������������, �� ������ �� ��������
        int iRelN = (index + i)%length;
        if (edgeList[edges[iRelN]].color == yellow) {
            index = iRelN;
            break;
        }
    }
    return edges[index];
}

int WF(int iE1, int iV, int iE2) {// iE1, iE2 - ������� ����� ������������� � ������� � �������� iV
    // fout << "WF: iE1= " << iE1 << " iV= " << iV << " iE2= " << iE2 << endl;
    vertexList[iV].color = black;
    int iVn = iV;
    while (iE1 != iE2) {
        edgeList[iE2].color = yellow;
        iVn = otherVertex(iVn, iE2);
        if (vertexList[iVn].color == green) {
            int iVnn = iVn;
            while (iVn != iV) {
                edgeList[iE2].color = green;
                iVn =otherVertex(iVn, iE2);
                Vertex &vertexiVn = vertexList[iVn];// ���� � �������� �����������
                vertexiVn.color = green;
                vector <int> &edges = vertexiVn.edges;
                iE2 = NextYellowEdge(edges,IndexOfEdge(edges,iE2));
                //fout << "WF: ww  iE2= " << iE2  << endl;
            }
            return iVnn;
        }
        vector <int> &iVnEdges=vertexList[iVn].edges;
        iE2 = NextEdge(iVnEdges, IndexOfEdge(iVnEdges,iE2));
        //fout << "WF: w  iE2= " << iE2 << endl;
        //printvertex(iVn);
    }
    vertexList[iV].color = green;
    return iV;
}

int IndexOfPrevGreenEdge(vector<int> &edges, int index){// ������� ���������� ������ ������� �������� ����� ������ ������� �������
    //����� ����� � �������� index � ������ ����� edges, ���� ������ �������� ����� ���, �� ������� �� ������ ������
    for (size_t i=1;i<edges.size();i++)
    {
        int i2=(index+i)%edges.size();
        if (edgeList[edges[i2]].color == green) {
            index = i2;
            break;
        }
    }
    return index;
}

int deep = 0;
int maxdeep = 0;
int countng = 0;

void NG(int iE, int iV1, int iV2)//iE-������� �����, iV1- ����������� ��� ������� �������, iV2- ������� �������, � ������� ��������� ��������,
//iV2 ����� ��������� � iV1; ���� �� ���������, �� ��������������, ��� ���������� ���� �� ������� ����� �� iV1 �� iV2. ����� ����� ����������� ��������
{
    if (iV1 == iV2) return;
    vector <int> &iV1Edges=vertexList[iV1].edges;
    int indiE = IndexOfEdge(iV1Edges,iE);
    int indiEn = IndexOfPrevGreenEdge(iV1Edges,indiE);
    // �� �������� ���������� ��������� �����������, ��� ����� � �������� iE ����� ����� ��������� ������� ������
    int iEn = iV1Edges[indiEn];
    int iV1n = otherVertex(iV1, iEn);
    int iEnn = NextEdge(iV1Edges, indiEn);
    if (iEnn == iE) NG(iEn, iV1n, iV2);//���� ������ �� ������� ������� �.8
    else {
        int iV2n = WF(iE, iV1, iEnn); //E'' ��� ������ ����� ����� E' � E ����� ����������� � ������� ����, ���� ��� ����������� ������� ����, ����������� ������� V1 � ������ ������� ��������
        if (iV2n == iV1) NG(iEn, iV1n, iV2);//��������� ������� ����, ������� ����� ��� ������ ������� ������ � ������� ������. �.9.2
        else {// ����� ����� ����� E � E' ���������� ������� ����, ������������� ������ � 9.3
            NG(iEn, iV1n, iV2n);//���� ������ �� ������� �������, �� �� ������� V'2, ������� ������� WF �9.3.2
            NG(iE, iV1, iV2);//����� NG � ���� �� ����������� ������������ �����, ������� ����� E � ��������� ������� ������� ����� �. 9.3.1
        }
    }
}

void AddEdgeWithVertices(int v1, int v2) {// ����������� � ������ edgeList ������ ������ (v1,v2)
    //����� ����� ����� ����������� � ������ ����������� ����� ��� ������ v1 � v2
    vertexList[v1].edges.push_back(edgeList.size());
    vertexList[v2].edges.push_back(edgeList.size());
    edgeList.push_back(Edge(v1, v2, black));
}   //-----------------------------------------------------------
#ifdef DEBUG
void printgrid(vector<vector<int>>&grid) {
    for (size_t i = 0; i < grid.size(); i++) {
        for (size_t j = 0; j < grid[i].size(); j++) {
            fout << grid[i][j];
            if (j < grid[i].size() - 1)  fout << " ";
        }
        fout << "\n";
    }
}

void printvertexList() {
    for (size_t i = 0; i < vertexList.size(); i++)
        printvertex(i);
}

void printedge(int i) {
    Edge& E = edgeList[i];
    fout << i << " " << E.v1 << " " << E.v2 << " " << E.color<< endl;
}

void printedgeList() {
    for (size_t i = 0; i < edgeList.size(); i++)
        printedge(i);
}
#endif // DEBUG

int main(int argc,char**argv){
    int64_t t1, t2;
    int flagout=1;
    vector<vector<int>>grid;// ��������� ������ ��� ������ ������� �� �����
    if(argc>1){
        inputFileName=argv[1];
        if(argc>2)
            resultFileName=argv[2];
        if(argc>3)
            sscanf(argv[3],"%d",&flagout);
    }
    ifstream fi(inputFileName);
    string s;
//������ �������� �����
    for (int i=0;getline(fi,s);i++)
    {
        stringstream ss;
        grid.push_back(vector<int>(0));
        ss << s;
        int tmp;
        while (ss >> tmp)
            grid[i].push_back(tmp);
    }
//���������� ������ � ������ ����������� �������
    int grid_width = grid[0].size();
	t1=GetTicks();
    // ��������� ����
    /*vector<int> buss0(grid_width, 0);
    vector<int> buss1(grid_width, 1);
    buss0[0] = 1;// ���������� ����������� ���� � ������������ ����
    grid.insert(grid.begin(), buss1);
    grid.insert(grid.begin(), buss0);
    grid.push_back(buss1);
    grid.push_back(buss0);*/

    int Lver = grid.size();
    int Lhor = grid[0].size();
    vector<vector<int>> grid_help; // ��� ����� ��������� ������� ������ � vertexList
//������������ vertexList � ���������������� ���������� ������� grid_help � ��������� ������ � vertexList
    for (int i = 0; i < Lver; i++) {
        grid_help.push_back(vector <int>(0));
        for (int j = 0; j < Lhor; j++)
            if (grid[i][j] > 0) {
                grid_help[i].push_back(vertexList.size());
                vertexList.push_back(Vertex(i, j, black));
            }
            else
                grid_help[i].push_back(grid[i][j]);
    }
// ������������ ������ ���� ����� � ���������� � ������ ������� ������ ����������� �� �����
// ������� ���������� ����������� ����� ������������� ��������� ������ ������� ������� ��� ������ �������
    for (size_t v1 = 0; v1 < vertexList.size(); v1++) {
        Vertex &vertex = vertexList[v1];
        int row = vertex.row;
        int col = vertex.col;
        if (row < Lver - 1) {
            size_t v2 = grid_help[row + 1][col];
            if (v2 > v1&& grid[row + 1][col] > 0) AddEdgeWithVertices(v1, v2);
        }
        if (col < Lhor - 1) {
            size_t v2 = grid_help[row][col + 1];
            if (v2 > v1&& grid[row][col + 1] > 0) AddEdgeWithVertices(v1, v2);
        }
    }
    //��������� Vin � Vout � ����� ������
    int Vin = vertexList.size();
    vertexList.push_back(Vertex(-2, -1, green));
    int Vout =vertexList.size();
    vertexList.push_back(Vertex(-1, -1, green));
    // � ��������� �������� ����������� ��� ������� ����� (Vout, Vin), (Vin, Vout)
    // ���������� ����� ��������� ���������� ������ ������� ������� ��� ������ �������
    int aE1 = edgeList.size(); // ������ ������ (�������) ����� (Vout, Vin) � edgeList
    edgeList.push_back(Edge(Vout, Vin, green));// ������� ������� �����
    vertexList[Vin].edges.push_back(aE1);
    vertexList[Vout].edges.push_back(aE1);
    for (int j = 0; j < Lhor; j++)//���� ��� ���������� ����� (�������) ����� �� Vout
        if (grid[0][j] > 0)
            AddEdgeWithVertices(Vout, grid_help[0][j]);
    int aE2 = edgeList.size(); // ������ ������� (��������) ����� (Vin, Vout) � edgeList (����� E0 � ����������)
    edgeList.push_back(Edge(Vin, Vout, green));// ������� ������� ������
    vertexList[Vin].edges.push_back(aE2);
    vertexList[Vout].edges.push_back(aE2);
    for (int j = Lhor-1; j >=0 ; j--)//���� ��� ���������� ������ (������) ����� �� Vin
        if (grid[Lver-1][j] > 0){
                AddEdgeWithVertices(Vin, grid_help[Lver - 1][j]);
                if(j<Lhor-1 && grid[Lver-1][j+1])//���� ���� ����� ������ �� �������, �� �������� ��������� ��� �����
                {
                    vector <int> &edges = vertexList[grid_help[Lver - 1][j]].edges;
                    iter_swap (edges.end()-1,edges.end()-2);
                }
            }

    /*printgrid(grid_help);
    fout<<"vertexList: "<<endl;
    printvertexList();
    fout<<"edgeList: "<<endl;
    printedgeList();*/
    NG(aE2, Vin, Vout);

	int bbmas=0;
	for (int i = 0; i < Lver; i++) {
        for (int j = 0; j < Lhor; j++) {
            if (grid[i][j] > 0) {
                int vI = grid_help[i][j];
                if (vertexList[vI].color != green) grid[i][j] = 0;
            }
            bbmas+=grid[i][j];
        }
    }

	t2=GetTicks();
	FILE *f3;
	char* repFile = (char*) "ReportAkhun_NoSort_Tail_NoBus.txt";
	f3=fopen(repFile,"a+");
	fseek(f3,0,SEEK_END);
	if(ftell(f3)==0)
		fprintf(f3,"Filename Lver Lhor Black_c Backbone_c Time \n");
	int q_backbone=bbmas, q_black= Vin;
	double sec = double (t2 - t1)/GetTicksPerSec();
	fprintf(f3,"%s %d %d %6d %6d %.7f\n",inputFileName,Lver,Lhor,q_black,q_backbone,sec);
	fclose(f3);

    printf("backbone mass= %d\n", bbmas);
    printf("runtime = %.7f\n",sec);
    if(flagout)
    {
        ofstream fo (resultFileName);
        for (int i = 0; i < Lver; i++) {
            for (int j = 0; j < Lhor; j++) {
                fo << grid[i][j];
                if (j < grid_width - 1)  fo << " ";
            }
            fo << "\n";
        }
    }

    return 0;
}
