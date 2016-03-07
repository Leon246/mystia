package Console;

/*step 3
 * 
 */

import java.io.*;

public class Sequence_muti_distribution {

	/**
	 * @param args
	 */
	public long[] stat;
	public static int dim = 2;
	
	public int coding(char aa) {
		int value = 0;
		switch (aa) {
        case 'G':
        	value = 0;
            break;
        case 'A':
        	value = 1;
            break;
        case 'P':
        	value = 2;
            break;
        case 'V':
        	value = 3;
            break;
        case 'L':
        	value = 4;
            break;
        case 'I':
        	value = 5;
            break;
        case 'M':
        	value = 6;
            break;
        case 'F':
        	value = 7;
            break;
        case 'Y':
        	value = 8;
            break;
        case 'W':
        	value = 9;
            break;
        case 'S':
        	value = 10;
            break;
        case 'T':
        	value = 11;
            break;
        case 'C':
        	value = 12;
            break;
        case 'N':
        	value = 13;
            break;
        case 'Q':
        	value = 14;
            break;
        case 'K':
        	value = 15;
            break;
        case 'H':
        	value = 16;
             break;
        case 'R':
        	value = 17;
            break;
        case 'D':
        	value = 18;
            break;
        case 'E':
        	value = 19;
            break;
        default:
        	System.out.print("aa coding error! \n");
        	break;
        }
		return value;
	}
	public char anti_coding(int code) {
		char aa = 'Z';
		switch(code) {
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
			System.out.print("aa anti_coding error \n");
		}
		return aa;
	}
	
	public Sequence_muti_distribution() throws IOException {
		String infilename =
//        "D:\\My Documents\\GThesis\\DatabaseBAK\\std_aa.xml";
          "D:\\My Documents\\GThesis\\DatabaseBAK\\std_aa_1M_01.xml";
//        "E:\\My Documents\\GThesis\\DatabaseBAK\\std_aa.xml";
//        "E:\\My Documents\\GThesis\\DatabaseBAK\\std_aa_1M_01.xml";

	    byte data[] = new byte[1];
	    char line[] = new char[120];
	    int i = 0, j = 0;
	    int k = 0;
	    int dim_para = 0;
	    int position = 0;
	    
	    stat = new long[(int)Math.pow(20, dim)];
	    for (k = 0; k < (int)Math.pow(20, dim); k ++) {
	    	stat[k] = 0;
	    }
		
	    try {
	        FileInputStream infile = new FileInputStream(infilename);
	        System.out.print("Reading data...\n");
	
	        while ((infile.read(data, 0, 1)) != -1) {
	            if (i > 14 && i < 135 && data[0] > 64 && data[0] < 90) {
	                line[j] = (char) data[0];
	                j++;
	            } else if (i == 138 && data[0] == 13 && j == 120) {
	            	for (k = 0; k < 121 - dim; k ++) {
	            		while (dim_para < dim) {
	            		//sequence 'ABC' => A * 20 ^ 0 + B * 20 ^ 1 + C * 20 ^ 2
	            			position += 
	            				coding(line[k + dim_para]) * (int)(Math.pow(20, dim_para));
	            			dim_para ++;
	            		}
	            		stat[position] ++;
	            		position = 0;
	            		dim_para = 0;
	            	}
	            }
	            i++;
	            if (data[0] == 10) {
	            	i = 0;
	            	j = 0;
	            }
	        }
	        System_out();
	        infile.close();
	    } catch (IOException e) {
	       throw e;
	    }
	}
	
	public void System_out() throws IOException {
		String outfilename = "D:\\My Documents\\GThesis\\Res\\sequence.txt";
		String str_filehead = ">|std_aa stat. with continues aa of : ";
		String strbuffer = String.valueOf(dim);
		str_filehead += strbuffer;
		
		if (!(dim > 0)) {
			System.out.print("Fatal error, please check dim set! \n");
			System.exit(1);
		}
		
		try {
	        FileOutputStream outfile = new FileOutputStream(outfilename);
	        DataOutputStream outdata = new DataOutputStream(outfile);
	        
	        outdata.writeUTF(str_filehead);
			outfile.write(32);
			outfile.write(13);
			int k = 0;
			int k_para = k;
			int dim_para = 0;
			int code[] = new int[dim];
			int MAX = (int)Math.pow(20, dim);
			
			//attention: in_data can not be overflowed
			//attention: in_data current define : 10
			final int in_data_dim = 10;
			byte in_data[] = new byte[in_data_dim];
			long value = 0;
			int loop_input = 0;
			
			for (k = 0; k < MAX; k ++) {
				outfile.write(10);
				k_para = k;
				for (dim_para = 0; dim_para < dim; dim_para ++) {
					code[dim - dim_para - 1] = k_para / (int)Math.pow(20, (dim - dim_para - 1));
					k_para %= (int)Math.pow(20,(dim - dim_para - 1));
				}
				for (dim_para = 0; dim_para < dim; dim_para ++) {
					outfile.write((byte)(anti_coding(code[dim_para])));
				}
				outfile.write(32);
				outfile.write(58);
				outfile.write(32);

				value = stat[k];
				strbuffer = (new Long(value)).toString();
				for (loop_input = 0; loop_input < strbuffer.length(); loop_input ++)
					outfile.write((byte)strbuffer.charAt(loop_input));
				
				outfile.write(in_data);
				outfile.write(13);
			}
			outdata.close();
			outfile.close();
			
			System.out.print("Finished\n");
		}catch(IOException e) {
			throw e;
		}
	}
	
	public static void main(String[] arg) throws IOException {
		new Sequence_muti_distribution();
	}
}

