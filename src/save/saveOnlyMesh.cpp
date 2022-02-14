#include <fstream>
#include <cstring>
#include <mpi.h>
using namespace std;

#include "save.h" 

void SEMO_Mesh_Save::vtu(SEMO_Mesh_Builder &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (.vtu format) ... ";
	}
	ofstream outputFile;
	string filenamePlot = "./save/plot." + to_string(rank) + ".vtu";
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.listPoints.size() << "\" NumberOfCells=\"" << mesh.listCells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	outputFile << "    </PointData>" << endl;
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	// Points
	outputFile << "    <Points>" << endl;
	// }
	outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	stringstream streamXYZ;
	// for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
	for(auto& point : mesh.points){
		outputFile << scientific << point.x << " " << point.y << " " << point.z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	for(auto& cell : mesh.cells){
		for(auto i : cell.points){
			outputFile << i << " ";
		}
		outputFile << endl;
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// offsets (cell's points offset)
	int cellFaceOffset = 0;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	cellFaceOffset = 0;
	for(auto& cell : mesh.cells){
		cellFaceOffset += cell.points.size();
		outputFile << cellFaceOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// types (cell's type, 42 = polyhedron)
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	for(auto& cell : mesh.cells){
		outputFile << "42" << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// outputFile << mesh.faces.size() << endl;
	for(auto& cell : mesh.cells){
		outputFile << cell.faces.size() << endl;
		for(auto& i : cell.faces){
			outputFile << mesh.faces[i].points.size() << " ";
			for(auto& j : mesh.faces[i].points){
				outputFile << j << " ";
			}
			outputFile << endl;
		}
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// faceoffsets (cell's face offset)
	int cellFacePointOffset = 0;
	
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	
	cellFacePointOffset = 0;
	for(auto& cell : mesh.cells){
		int numbering = 1 + cell.faces.size();
		for(auto& i : cell.faces){
			numbering += mesh.faces[i].points.size();
		}
		cellFacePointOffset += numbering;
		outputFile << cellFacePointOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	

	// additional informations
	outputFile << " <owner>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.owner << " ";
	}
	outputFile << endl;
	outputFile << " </owner>" << endl;
	
	outputFile << " <neighbour>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.neighbour << " ";
	}
	outputFile << endl;
	outputFile << " </neighbour>" << endl;
	
	outputFile << " <bcName>" << endl;
	for(auto& boundary : mesh.boundary){
		// cout << boundary.name << endl;
		outputFile << boundary.name << " ";
	}
	outputFile << endl;
	outputFile << " </bcName>" << endl;
	
	outputFile << " <bcStartFace>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.startFace << " ";
	}
	outputFile << endl;
	outputFile << " </bcStartFace>" << endl;
	
	outputFile << " <bcNFaces>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.nFaces << " ";
	}
	outputFile << endl;
	outputFile << " </bcNFaces>" << endl;
	
	outputFile << " <bcNeighbProcNo>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.neighbProcNo << " ";
	}
	outputFile << endl;
	outputFile << " </bcNeighbProcNo>" << endl;
	
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	// }
	
	
	
}




void SEMO_Mesh_Save::vtu(string folder, int myRank, SEMO_Mesh_Builder &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	
	char folder_name[1000];
	strcpy(folder_name, folder.c_str());
	mkdirs(folder_name);
	
	// fs::path Paths(folder);
	// if( ! exists(Paths) ){
		// create_directories(Paths);
	// }
	
	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(myRank) + ".vtu";
	

	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << filenamePlot << ") ... ";
	}
	
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.listPoints.size() << "\" NumberOfCells=\"" << mesh.listCells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	outputFile << "    </PointData>" << endl;
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	// Points
	outputFile << "    <Points>" << endl;
	// }
	outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	stringstream streamXYZ;  
	outputFile.precision(20);  
	// for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
	for(auto& point : mesh.points){
		outputFile << scientific << point.x << " " << point.y << " " << point.z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	for(auto& cell : mesh.cells){
		for(auto i : cell.points){
			outputFile << i << " ";
		}
		outputFile << endl;
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// offsets (cell's points offset)
	int cellFaceOffset = 0;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	cellFaceOffset = 0;
	for(auto& cell : mesh.cells){
		cellFaceOffset += cell.points.size();
		outputFile << cellFaceOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// types (cell's type, 42 = polyhedron)
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	for(auto& cell : mesh.cells){
		outputFile << "42" << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// outputFile << mesh.faces.size() << endl;
	for(auto& cell : mesh.cells){
		outputFile << cell.faces.size() << endl;
		for(auto& i : cell.faces){
			outputFile << mesh.faces[i].points.size() << " ";
			for(auto& j : mesh.faces[i].points){
				outputFile << j << " ";
			}
			outputFile << endl;
		}
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// faceoffsets (cell's face offset)
	int cellFacePointOffset = 0;
	
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	
	cellFacePointOffset = 0;
	for(auto& cell : mesh.cells){
		int numbering = 1 + cell.faces.size();
		for(auto& i : cell.faces){
			numbering += mesh.faces[i].points.size();
		}
		cellFacePointOffset += numbering;
		outputFile << cellFacePointOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	

	// additional informations
	outputFile << " <owner>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.owner << " ";
	}
	outputFile << endl;
	outputFile << " </owner>" << endl;
	
	outputFile << " <neighbour>" << endl;
	for(auto& face : mesh.faces){
		outputFile << face.neighbour << " ";
	}
	outputFile << endl;
	outputFile << " </neighbour>" << endl;
	
	outputFile << " <bcName>" << endl;
	for(auto& boundary : mesh.boundary){
		// cout << boundary.name << endl;
		outputFile << boundary.name << " ";
	}
	outputFile << endl;
	outputFile << " </bcName>" << endl;
	
	outputFile << " <bcStartFace>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.startFace << " ";
	}
	outputFile << endl;
	outputFile << " </bcStartFace>" << endl;
	
	outputFile << " <bcNFaces>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.nFaces << " ";
	}
	outputFile << endl;
	outputFile << " </bcNFaces>" << endl;
	
	outputFile << " <bcNeighbProcNo>" << endl;
	for(auto& boundary : mesh.boundary){
		outputFile << boundary.neighbProcNo << " ";
	}
	outputFile << endl;
	outputFile << " </bcNeighbProcNo>" << endl;
	
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	// }
	
	
	
}














































// void MASCH_Mesh_Save::vtu(MASCH_Mesh &mesh){
	
	// int rank = MPI::COMM_WORLD.Get_rank();
	
	// // if(rank==0){
		// // cout << "┌────────────────────────────────────────────────────" << endl;
		// // cout << "| execute save file (.vtu format) ... ";
	// // }
	// // ofstream outputFile;
	// // string filenamePlot = "./save/plot." + to_string(rank) + ".vtu";
	// // outputFile.open(filenamePlot);
	// // if(outputFile.fail()){
		// // cerr << "Unable to write file for writing." << endl;
		// // MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	// // }
	
	// // // string out_line;
	// // outputFile << "<?xml version=\"1.0\"?>" << endl;
	// // outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	// // outputFile << "  <UnstructuredGrid>" << endl;
	// // outputFile << "   <Piece NumberOfPoints=\"" << mesh.listPoints.size() << "\" NumberOfCells=\"" << mesh.listCells.size() << "\">" << endl;
	
	// // // Points data
	// // outputFile << "    <PointData>" << endl;
	// // outputFile << "    </PointData>" << endl;
	// // // Cells data
	// // outputFile << "    <CellData>" << endl;
	// // outputFile << "    </CellData>" << endl;
	// // // Points
	// // outputFile << "    <Points>" << endl;
	// // // }
	// // outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	// // stringstream streamXYZ;
	// // // for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
	// // for(auto& point : mesh.points){
		// // outputFile << scientific << point.x << " " << point.y << " " << point.z << endl;

	// // }
	
	// // outputFile << "    </DataArray>" << endl;
	// // outputFile << "   </Points>" << endl;
	
	// // // cells
	// // outputFile << "   <Cells>" << endl; 
	// // // connectivity (cell's points)
	// // outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	// // for(auto& cell : mesh.cells){
		// // for(auto i : cell.points){
			// // outputFile << i << " ";
		// // }
		// // outputFile << endl;
	// // }
	
	// // outputFile << "    </DataArray>" << endl;
	
	// // // offsets (cell's points offset)
	// // int cellFaceOffset = 0;
	// // outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	// // cellFaceOffset = 0;
	// // for(auto& cell : mesh.cells){
		// // cellFaceOffset += cell.points.size();
		// // outputFile << cellFaceOffset << " ";
	// // }
	// // outputFile << endl;
	
	// // outputFile << "    </DataArray>" << endl;
	
	// // // types (cell's type, 42 = polyhedron)
	// // outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	// // for(auto& cell : mesh.cells){
		// // outputFile << "42" << " ";
	// // }
	// // outputFile << endl;
	
	// // outputFile << "    </DataArray>" << endl;
	
	// // // faces (cell's faces number, each face's point number, cell's faces's points)
	// // outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// // // outputFile << mesh.faces.size() << endl;
	// // for(auto& cell : mesh.cells){
		// // outputFile << cell.faces.size() << endl;
		// // for(auto& i : cell.faces){
			// // outputFile << mesh.faces[i].points.size() << " ";
			// // for(auto& j : mesh.faces[i].points){
				// // outputFile << j << " ";
			// // }
			// // outputFile << endl;
		// // }
	// // }
	
	// // outputFile << "    </DataArray>" << endl;
	
	// // // faceoffsets (cell's face offset)
	// // int cellFacePointOffset = 0;
	
	// // outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	
	// // cellFacePointOffset = 0;
	// // for(auto& cell : mesh.cells){
		// // int numbering = 1 + cell.faces.size();
		// // for(auto& i : cell.faces){
			// // numbering += mesh.faces[i].points.size();
		// // }
		// // cellFacePointOffset += numbering;
		// // outputFile << cellFacePointOffset << " ";
	// // }
	// // outputFile << endl;
	
	// // outputFile << "    </DataArray>" << endl;
	// // outputFile << "   </Cells>" << endl;
	
	
	// // outputFile << "  </Piece>" << endl;
	// // outputFile << " </UnstructuredGrid>" << endl;
	

	// // // additional informations
	// // outputFile << " <owner>" << endl;
	// // for(auto& face : mesh.faces){
		// // outputFile << face.owner << " ";
	// // }
	// // outputFile << endl;
	// // outputFile << " </owner>" << endl;
	
	// // outputFile << " <neighbour>" << endl;
	// // for(auto& face : mesh.faces){
		// // outputFile << face.neighbour << " ";
	// // }
	// // outputFile << endl;
	// // outputFile << " </neighbour>" << endl;
	
	// // outputFile << " <bcName>" << endl;
	// // for(auto& boundary : mesh.boundary){
		// // // cout << boundary.name << endl;
		// // outputFile << boundary.name << " ";
	// // }
	// // outputFile << endl;
	// // outputFile << " </bcName>" << endl;
	
	// // outputFile << " <bcStartFace>" << endl;
	// // for(auto& boundary : mesh.boundary){
		// // outputFile << boundary.startFace << " ";
	// // }
	// // outputFile << endl;
	// // outputFile << " </bcStartFace>" << endl;
	
	// // outputFile << " <bcNFaces>" << endl;
	// // for(auto& boundary : mesh.boundary){
		// // outputFile << boundary.nFaces << " ";
	// // }
	// // outputFile << endl;
	// // outputFile << " </bcNFaces>" << endl;
	
	// // outputFile << " <bcNeighbProcNo>" << endl;
	// // for(auto& boundary : mesh.boundary){
		// // outputFile << boundary.neighbProcNo << " ";
	// // }
	// // outputFile << endl;
	// // outputFile << " </bcNeighbProcNo>" << endl;
	
	
	// // outputFile << "</VTKFile>" << endl;
	
	// // outputFile.close();
	// // if(rank==0){
		// // cout << "-> completed" << endl;
		// // cout << "└────────────────────────────────────────────────────" << endl;
	// // }
	
	// // // }
	
	
	
// }




void MASCH_Mesh_Save::vtu(string folder, int myRank, MASCH_Mesh &mesh){
	
	int rank = MPI::COMM_WORLD.Get_rank();
	
	// 폴더 만들기
    auto ret = filesystem::create_directories(folder);
	
	ofstream outputFile;
	string filenamePlot = folder + "plot." + to_string(myRank) + ".vtu";
	
	if(rank==0){
		cout << "┌────────────────────────────────────────────────────" << endl;
		cout << "| execute save file (" << filenamePlot << ") ... ";
	}
	
	outputFile.open(filenamePlot);
	if(outputFile.fail()){
		cerr << "Unable to write file for writing." << endl;
		MPI_Abort(MPI_COMM_WORLD,EXIT_FAILURE);
	}
	
	// string out_line;
	outputFile << "<?xml version=\"1.0\"?>" << endl;
	outputFile << " <VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">" << endl;
	outputFile << "  <UnstructuredGrid>" << endl;
	outputFile << "   <Piece NumberOfPoints=\"" << mesh.points.size() << "\" NumberOfCells=\"" << mesh.cells.size() << "\">" << endl;
	
	// Points data
	outputFile << "    <PointData>" << endl;
	outputFile << "    </PointData>" << endl;
	// Cells data
	outputFile << "    <CellData>" << endl;
	outputFile << "    </CellData>" << endl;
	// Points
	outputFile << "    <Points>" << endl;
	// }
	outputFile << "     <DataArray type=\"Float64\" Name=\"NodeCoordinates\" NumberOfComponents=\"3\" format=\"ascii\">" << endl;

	stringstream streamXYZ;  
	outputFile.precision(20);  
	// for(auto iter=mesh.points.begin(); iter!=mesh.points.end(); iter++){
	for(auto& point : mesh.points){
		outputFile << scientific << point.x << " " << point.y << " " << point.z << endl;

	}
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Points>" << endl;
	
	// cells
	outputFile << "   <Cells>" << endl; 
	// connectivity (cell's points)
	outputFile << "    <DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">" << endl;

	for(auto& cell : mesh.cells){
		for(auto i : cell.ipoints){
			outputFile << i << " ";
		}
		outputFile << endl;
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// offsets (cell's points offset)
	int cellFaceOffset = 0;
	outputFile << "    <DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">" << endl;
	
	cellFaceOffset = 0;
	for(auto& cell : mesh.cells){
		cellFaceOffset += cell.ipoints.size();
		outputFile << cellFaceOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// types (cell's type, 42 = polyhedron)
	outputFile << "    <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">" << endl;
	
	for(auto& cell : mesh.cells){
		outputFile << "42" << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	
	// faces (cell's faces number, each face's point number, cell's faces's points)
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faces\" format=\"ascii\">" << endl;
	
	// outputFile << mesh.faces.size() << endl;
	for(auto& cell : mesh.cells){
		outputFile << cell.ifaces.size() << endl;
		for(auto& i : cell.ifaces){
			outputFile << mesh.faces[i].ipoints.size() << " ";
			for(auto& j : mesh.faces[i].ipoints){
				outputFile << j << " ";
			}
			outputFile << endl;
		}
	}
	
	outputFile << "    </DataArray>" << endl;
	
	// faceoffsets (cell's face offset)
	int cellFacePointOffset = 0;
	
	outputFile << "    <DataArray type=\"Int64\" IdType=\"1\" Name=\"faceoffsets\" format=\"ascii\">" << endl;
	
	cellFacePointOffset = 0;
	for(auto& cell : mesh.cells){
		int numbering = 1 + cell.ifaces.size();
		for(auto& i : cell.ifaces){
			numbering += mesh.faces[i].ipoints.size();
		}
		cellFacePointOffset += numbering;
		outputFile << cellFacePointOffset << " ";
	}
	outputFile << endl;
	
	outputFile << "    </DataArray>" << endl;
	outputFile << "   </Cells>" << endl;
	
	
	outputFile << "  </Piece>" << endl;
	outputFile << " </UnstructuredGrid>" << endl;
	
	
	
	

	// additional informations
	{
		string saveFormat = "ascii";
		outputFile << " <DataArray type=\"Int32\" Name=\"owner\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		for(auto& face : mesh.faces){
			// values.push_back(face.owner);
			outputFile << face.iL << " ";
		}
		// writeDatasAtVTU(controls, outputFile, values);
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
	}
	{
		string saveFormat = "ascii";
		outputFile << " <DataArray type=\"Int32\" Name=\"neighbour\" format=\"" << saveFormat << "\">" << endl;
		// vector<int> values;
		for(auto& face : mesh.faces){
			// values.push_back(face.neighbour);
			if(face.getType() == MASCH_Face_Types::INTERNAL){
				outputFile << face.iR << " ";
			}
		}
		// writeDatasAtVTU(controls, outputFile, values);
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
	}
	
	

	// boundary informations
	{
		outputFile << " <DataArray type=\"Char\" Name=\"bcName\" format=\"" << "ascii" << "\">" << endl;
		// for(auto& boundary : mesh.boundary){
			// // cout << boundary.name << endl;
			// // trim;
			// string bcName = boundary.name;
			
			// bcName.erase(std::find_if(bcName.rbegin(), bcName.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), bcName.end());
			// bcName.erase(bcName.begin(), std::find_if(bcName.begin(), bcName.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
			
			
			// // outputFile << boundary.name << " ";
			// outputFile << bcName << " ";
		// }
		for(auto& boundary : mesh.boundaries){
			// cout << boundary.name << endl;
			outputFile << boundary.name << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcStartFace\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.startFace << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNFaces\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			outputFile << boundary.nFaces << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"bcNeighbProcNo\" format=\"" << "ascii" << "\">" << endl;
		for(auto& boundary : mesh.boundaries){
			// cout << boundary.rightProcNo << endl;
			outputFile << boundary.rightProcNo << " ";
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		outputFile << " <DataArray type=\"Int32\" Name=\"connPoints\" format=\"" << "ascii" << "\">" << endl;
		for(int i=0; i<mesh.points.size(); ++i){
			auto& point = mesh.points[i];
			for(auto& [item, item2] : point.connPoints){
				outputFile << i << " " << item << " " << item2 << " ";
			}
			// if(!point.connPoints.empty()) outputFile << endl;
		}
		outputFile << endl;
		outputFile << " </DataArray>" << endl;
		
		
	}
	
	outputFile << "</VTKFile>" << endl;
	
	outputFile.close();
	if(rank==0){
		cout << "-> completed" << endl;
		cout << "└────────────────────────────────────────────────────" << endl;
	}
	
	// }
	
	
	
}


