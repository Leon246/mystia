package Console;

import java.io.*;

public class Sequence_std_out {

	/**
	 * @param args
	 */
	
	public static int dim = 3;
	
	public char anti_coding(int code) {
		char aa = 'Z';
		switch (code){
		case 0:
			aa = 'G';
			break;
		case 1:
			aa = 'A';
			break;
		case 2:
			aa = 'P';
			break;
		case 3:
			aa = 'V';
			break;
		case 4:
			aa = 'L';
			break;
		case 5:
			aa = 'I';
			break;
		case 6:
			aa = 'M';
			break;
		case 7:
			aa = 'F';
			break;
		case 8:
			aa = 'Y';
			break;
		case 9:
			aa = 'W';
			break;
		case 10:
			aa = 'S';
			break;
		case 11:
			aa = 'T';
			break;
		case 12:
			aa = 'C';
			break;
		case 13:
			aa = 'N';
			break;
		case 14:
			aa = 'Q';
			break;
		case 15:
			aa = 'K';
			break;
		case 16:
			aa = 'H';
			break;
		case 17:
			aa = 'R';
			break;
		case 18:
			aa = 'D';
			break;
		case 19:
			aa = 'E';
			break;
		default:
			System.out.print("aa anti_coding error! \n");
			break;
		}
		return aa;
	}

	public Sequence_std_out() throws IOException{
		String outfilename = "D:\\My Documents\\GThesis\\Res\\sequence.txt";
		
		if (!(dim > 0)) {
			System.out.print("Fatal error, please check dim set! \n");
			System.exit(1);
		}
		
		try {
	        FileOutputStream outfile = new FileOutputStream(outfilename);
			
			int k = 0;
			int k_para = k;
			int dim_para = 0;
			int code[] = new int[dim];
			int MAX = (int)Math.pow(20, dim);
			
			for (k = 0; k < MAX; k ++) {
				k_para = k;
				for (dim_para = 0; dim_para < dim; dim_para ++) {
					code[dim - dim_para - 1] = k_para / (int)Math.pow(20, (dim - dim_para - 1));
					k_para %= (int)Math.pow(20,(dim - dim_para - 1));
				}
				for (dim_para = 0; dim_para < dim; dim_para ++) {
					outfile.write((byte)(anti_coding(code[dim_para])));
				}
				outfile.write(13);
				outfile.write(10);
			}
			outfile.write(3);
			outfile.close();
		}catch(IOException e) {
			throw e;
		}
	}
	
	public static void main(String[] args) throws IOException{
		// TODO Auto-generated method stub
		new Sequence_std_out();
	}

}