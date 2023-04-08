using System;
using System.Collections.Generic;
using System.IO;
using System.Text;
using MGroup.LinearAlgebra.Matrices;
using MGroup.MSolve.Discretization;
using MGroup.MSolve.Discretization.Entities;


namespace MGroup.DrugDeliveryModel.Tests.Integration
{
	public class VtkFileWriter2

	{
		public const string vtkReaderVersion = "4.1";
		private readonly StreamWriter writer;
		private bool writeFieldsNext;
		public int dimension;
		public int cellCode;
		//private int numVertices = -1;

		public VtkFileWriter2(string filePath, int dimension, CellType cellType)
		{
			this.dimension = dimension;
			if (cellType==CellType.Quad4)
			{
				cellCode = 9;
			}
			else if (cellType==CellType.Hexa8)
			{
				cellCode = 12;
			}
			else if (cellType==CellType.Tet4)
			{
				cellCode = 10;
			}
			this.writer = new StreamWriter(filePath);
			writer.Write("# vtk DataFile Version ");
			writer.WriteLine(vtkReaderVersion);
			//writer.WriteLine(filePath);
			writer.WriteLine("Header:");
			writer.Write("ASCII\n\n");
			writeFieldsNext = false;
		}

		public void Dispose()
		{
			if (writer != null) writer.Dispose();
		}

		//TODO: Perhaps the mesh should be injected into the contructor

		//Solid/void application
		public void WriteMesh(Dictionary<int, INode>  nodes, Dictionary<int, IElementType> elements)
		{
			if (writeFieldsNext) throw new InvalidOperationException("A mesh has already been written.");

			// Vertices 
			//this.numVertices = mesh.NumOutVertices;
			writer.WriteLine("DATASET UNSTRUCTURED_GRID");
			writer.WriteLine($"POINTS {nodes.Count} double");
			//if (nodes[0].Z != null)
			//{
			//	for (int i = 0; i < nodeCoords.NumRows; ++i)
			//	{
			//		writer.WriteLine($"{nodeCoords[i, 0]} {nodeCoords[i, 1]} {nodeCoords[i, 2]}");
			//	}
			//}
			//else
			//{
			if (dimension==2)
			{
				for (int i = 0; i < nodes.Count; ++i)
				{
					writer.WriteLine($"{nodes[i].X} {nodes[i].Y} 0.0");
				}
			}
			else if (dimension == 3)
			{
				for (int i = 0; i < nodes.Count; ++i)
				{
					writer.WriteLine($"{nodes[i].X} {nodes[i].Y} {nodes[i].Z}");
				}
			}
			//}

			// Cell connectivity
			int cellDataCount = 0;
			for (int i = 0; i < elements.Count; ++i)
			{
				cellDataCount += 1 + elements[i].Nodes.Count;
			}
			writer.WriteLine($"\nCELLS {elements.Count} {cellDataCount}");
			for (int i = 0; i < elements.Count; ++i)
			{
				writer.Write(elements[i].Nodes.Count);
				for (int j = 0; j < elements[i].Nodes.Count; ++j)
				{
					writer.Write(' ');
					writer.Write(elements[i].Nodes[j].ID);
				}
				writer.WriteLine();
			}

			// Cell types
			writer.WriteLine("\nCELL_TYPES " + elements.Count);
			for (int i = 0; i < elements.Count; ++i)
			{
				writer.WriteLine(cellCode);
			}
		}

		//Multiphase application
		//public void WriteMesh(Dictionary<int, Node> nodes, Dictionary<int, IElementType> elements)
		//{
		//	if (writeFieldsNext) throw new InvalidOperationException("A mesh has already been written.");

		//	// Vertices 
		//	//this.numVertices = mesh.NumOutVertices;
		//	writer.WriteLine("DATASET UNSTRUCTURED_GRID");
		//	writer.WriteLine($"POINTS {nodes.Count} double");
		//	//if (nodes[0].Z != null)
		//	//{
		//	//	for (int i = 0; i < nodeCoords.NumRows; ++i)
		//	//	{
		//	//		writer.WriteLine($"{nodeCoords[i, 0]} {nodeCoords[i, 1]} {nodeCoords[i, 2]}");
		//	//	}
		//	//}
		//	//else
		//	//{
		//	if (dimension == 2)
		//	{
		//		for (int i = 0; i < nodes.Count; ++i)
		//		{
		//			writer.WriteLine($"{nodes[i].X} {nodes[i].Y} 0.0");
		//		}
		//	}
		//	else if (dimension == 3)
		//	{
		//		for (int i = 0; i < nodes.Count; ++i)
		//		{
		//			writer.WriteLine($"{nodes[i].X} {nodes[i].Y} {nodes[i].Z}");
		//		}
		//	}
		//	//}

		//	// Cell connectivity
		//	int cellDataCount = 0;
		//	for (int i = 0; i < elements.Count; ++i)
		//	{
		//		cellDataCount += 1 + elements[i].Nodes.Count;
		//	}
		//	writer.WriteLine($"\nCELLS {elements.Count} {cellDataCount}");
		//	for (int i = 0; i < elements.Count; ++i)
		//	{
		//		writer.Write(elements[i].Nodes.Count);
		//		for (int j = 0; j < elements[i].Nodes.Count; ++j)
		//		{
		//			writer.Write(' ');
		//			writer.Write(elements[i].Nodes[j].ID);
		//		}
		//		writer.WriteLine();
		//	}

		//	// Cell types
		//	writer.WriteLine("\nCELL_TYPES " + elements.Count);
		//	for (int i = 0; i < elements.Count; ++i)
		//	{
		//		writer.WriteLine(cellCode);
		//	}
		//}


		public void WriteScalarField(string fieldName, IList<double> scalars)
		{
			WriteFieldsHeader(scalars.Count);
			writer.WriteLine($"SCALARS {fieldName} double 1");
			writer.WriteLine("LOOKUP_TABLE default");
			foreach (double value in scalars)
			{
				writer.WriteLine(value);
			};
			writer.WriteLine();
		}
		

		public void WriteVectorField(string fieldName, IList<double[]> vectors)
		{
			WriteFieldsHeader(vectors.Count);
			writer.WriteLine($"VECTORS {fieldName} double");
			foreach (double[] vector in vectors)
			{
				if (vector.Length == 3)
				{
					writer.WriteLine($"{vector[0]} {vector[1]} {vector[2]}");

				}
				else
				{
					writer.WriteLine($"{vector[0]} {vector[1]} 0.0");
				}
			}
			writer.WriteLine();
		}

		/// <summary>
		/// If the user only wants the mesh, this should not be called. Therefore only call it if one or more field output is 
		/// written.
		/// </summary>
		private void WriteFieldsHeader(int numVertices)
		{
			if (!writeFieldsNext) // Fields header
			{
				writer.Write("\n\n");
				writer.WriteLine("POINT_DATA " + numVertices);
				writeFieldsNext = true;
			}
		}
	}
}
