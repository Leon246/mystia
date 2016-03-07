package Console;

import java.io.*;
import java.util.Random;

public class DNA_artifical_generator {
	
	public int random_coding(double g, double c, double a, double t, Random r) {
		double na_value = r.nextDouble();
		c += g;
		a += c;
		t += a;
		if (na_value >= 0.0 && na_value < g) {
			return 71;
		}else if (na_value >= g && na_value < c) {
			return 67;
		}else if (na_value >= c && na_value < a) {
			return 65;
		}else if (na_value >= a) {
			return 84;
		}else {
			System.out.print("Internal Error\n");
			System.exit(1);
			return 0;
		}
	}
	
	public void original_fna(int file_num) throws IOException {
		String fna_filename = 
			"D:\\My Documents\\GThesis\\fna_data\\NC_";
		
		try {
			String temp = new String();
			temp = "0" + file_num;
			while (temp.length() < 6) {
				temp = "0" + temp;
			}
			
			fna_filename += temp + ".fna";
			
			FileOutputStream outfile = new FileOutputStream(fna_filename);
			DataOutputStream outdata = new DataOutputStream(outfile);
			
			String firstline = 
				">original random";
			String line = new String();
			outdata.writeBytes(firstline);
			outfile.write(13);
			outfile.write(10);

			
			int i = 0, j = 0;
			double value = 0;
			//value = 0.25 + (double)(file_num - 10) * 0.01;
			Random r = new Random();
			char insert;
			
			for (i = 0; i < 2000000; i ++) {
				
				insert = (char)random_coding(0.25, 0.25, 0.25, 0.25, r);
				line += insert;
				j ++;
				
				if (j == 70) {
					outfile.write(13);
					outdata.writeBytes(line);
					line = "";
					outfile.write(10);
					j = 0;
				}
				
				//not good
				
				if (insert == 'T') {
					insert = (char)random_coding(0.32, 0.32, 0.04, 0.32, r);
					line += insert;
					j ++;
				}
				if (j == 70) {
					outfile.write(13);					
					outdata.writeBytes(line);
					line = "";
					outfile.write(10);
					j = 0;
				}
				/*
				if (insert == 'T') {
					insert = (char)random_coding(0.3, 0.3, 0.1, 0.3, r);
					line += insert;
					j ++;
				}
				if (j == 70) {
					outfile.write(13);					
					outdata.writeBytes(line);
					line = "";
					outfile.write(10);
					j = 0;
				}

				if (insert == 'A') {
					insert = (char)random_coding(0.01, 0.33, 0.33, 0.33, r);
					line += insert;
					j ++;
				}
				if (j == 70) {
					outfile.write(13);					
					outdata.writeBytes(line);
					line = "";
					outfile.write(10);
					j = 0;
				}
				*/
			}
			outdata.close();
			outfile.close();
		}catch(IOException e) {
			throw e;
		}
	}
	
	public DNA_artifical_generator() throws IOException{
		int file_num = 41;
//		while (file_num < 19) {
			original_fna(file_num);
//			file_num ++;
//		}
	}

	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		new DNA_artifical_generator();
	}

}
