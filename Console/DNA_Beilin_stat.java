package Console;

/*Bei-lin Hao K_string
 *notice that each line begins with 10 (\r) and ends with 13(\n)
 *confirmed that no more error exits
 *read as loop
 *version 1.0
 */

import java.io.DataOutputStream;
import java.io.FileInputStream;
import java.io.DataInputStream;
import java.io.FileOutputStream;
import java.io.IOException;

public class DNA_Beilin_stat {

	public static int K_dim = 9;
	public static int start_point = 2;
	public long[] K_string;
	
	public int X_coding(char na) {
		int code = 2;
		switch (na) {
		case 'G':
            code = 0;
            break;
        case 'C':
            code = 1;
            break;
        case 'A':
            code = 0;
            break;
        case 'T':
            code = 1;
            break;
        default:
        	System.out.print("Nucleic acid X coding error \n");
        	break;
		}
		return code;
	}
	public int Y_coding(char na) {
		int code = 2;
		switch (na) {
		case 'G':
            code = 0;
            break;
        case 'C':
            code = 0;
            break;
        case 'A':
            code = 1;
            break;
        case 'T':
            code = 1;
            break;
		default:
			System.out.print("Nucleic acid Y coding error \n");
			break;
		}
		return code;
	}

	public char[] anti_coding(long K_value) {
		char k_string[] = new char[K_dim];
		char k_temp = 'U';
		int X = 0, Y = 0;
		int line = (int)Math.pow(2, K_dim);
		X = (int)(K_value % (long)line);
		Y = (int)(K_value / (long)line);
		
		int k_loop = 0;
		while (k_loop < K_dim) {
			line /= 2;
			if (X < line && Y < line) {
				k_temp = 'G';
			}else if (X >= line && Y < line) {
				k_temp = 'C';
				X -= line;
			}else if (X < line && Y >= line) {
				k_temp = 'A';
				Y -= line;
			}else if (X >= line && Y >= line) {
				k_temp = 'T';
				X -= line;
				Y -= line;
			}
			
			if (k_temp == 'U') {
				System.out.print("k_temp has not been used!\n");
				System.exit(2);
			}else {
				k_string[k_loop] = k_temp;
				k_loop ++;
			}
		}
		
		return k_string;
	}
	
	public DNA_Beilin_stat() throws IOException {

        String strfilename =
			"D:\\My Documents\\GThesis\\Data\\NC_000913.fna";
//          "D:\\My Documents\\GThesis\\Data\\Escherichia_coli_K12\\NC_000913.fna";
//    		"D:\\My Documents\\GThesis\\Data\\Bacillus_subtilis\\NC_000964.fna";
//    		"D:\\My Documents\\GThesis\\Data\\Streptococcus_mutans\\NC_004350.fna";
//    		"D:\\My Documents\\GThesis\\Data\\Thermoanaerobacter_tengcongensis\\NC_003869.fna";
//		    "E:\\My Documents\\GThesis\\Data\\NC_000913.fna";
//      "E:\\My Documents\\GThesis\\Data\\Escherichia_coli_K12\\NC_000913.fna";
//    	"E:\\My Documents\\GThesis\\Data\\Bacillus_subtilis\\NC_000964.fna";
//		"E:\\My Documents\\GThesis\\Data\\Streptococcus_mutans\\NC_004350.fna";
//		"E:\\My Documents\\GThesis\\Data\\Thermoanaerobacter_tengcongensis\\NC_003869.fna";

        //Define K_string, need not to be reset to zero
        int K = K_dim;
        int k = 0;
        int n = 0;
        K_string = new long[(int) Math.pow(4, (double)K)];

        try {
            FileInputStream infile = new FileInputStream(strfilename);
            DataInputStream indata = new DataInputStream(infile);

            int X = 0, Y = 0, value = 0;
            int half = (int)Math.pow(2., (double)(K - 1)), total = half * 2;
            byte data[] = new byte[1];
            char nacid;
            
            byte store[] = new byte[K_dim - 1];
            int store_loop = K_dim;
            
            boolean firstline = true;
            
            while ((indata.read(data, 0, 1)) != -1) {
            	//be very careful with the .fna format
	        	//jump over the first line only
	        	if (data[0] == 10 && firstline) {
	        		firstline = false;
	    	        System.out.print("First line passed! \n");
		        }else if (!firstline){
	                if ((data[0] == 10 || data[0] == 13) && !firstline) {
	                    continue;
	                } else {
	                	if (store_loop > 1) {
	                		store[K_dim - store_loop] = data[0];
	                		store_loop --;
	                	}
	                    //Total number of nucleic acid
	                    n++;
	                    //set nacid
	                    nacid = (char) data[0];
	                    if (k < K) {
	                        //X-axis of sequence
	                        value = X_coding(nacid);
	                        X += (int) Math.pow(2, K - k - 1) * value;
	                        //Y-axis of sequence
	                        value = Y_coding(nacid);
	                        Y += (int) Math.pow(2, K - k - 1) * value;
	
	                        k++;
	                        //Don't forget the first sequence
	                        if (k == K) {
	                            K_string[X + Y * total]++;
	                        }
	                    } else {
	                        //next nucleic acid
	                        value = X_coding(nacid);
	                        //next sequence
	                        if (X >= half) {
	                            X -= half;
	                        }
	                        X *= 2;
	                        X += value;
	                        value = Y_coding(nacid);
	                        //next sequence
	                        if (Y >= half) {
	                            Y -= half;
	                        }
	                        Y *= 2;
	                        Y += value;
	                        //save sequence
	                        K_string[X + Y * total]++;
	                    }
	                }
		        }else {
		        	continue;
		        }
            }
            //read as a loop
            for (store_loop = K_dim; store_loop > 1; store_loop --) {
            	nacid = (char)store[K_dim - store_loop];
            	//next nucleic acid
                value = X_coding(nacid);
                //next sequence
                if (X >= half) {
                    X -= half;
                }
                X *= 2;
                X += value;
                value = Y_coding(nacid);
                //next sequence
                if (Y >= half) {
                    Y -= half;
                }
                Y *= 2;
                Y += value;
                //save sequence
                K_string[X + Y * total]++;
            }
            
            System_out();
            //Total number of nucleic acid
            System.out.print("Total number : " + n + "\n");

            //close file
            indata.close();
            infile.close();
        } catch (IOException e) {
            throw e;
        }
	}
	
	public void System_out() throws IOException {
		if (K_dim == 3) {
			System.out.print("   TTT F Phe  " + K_string[63]);
	        System.out.print("   TCT S Ser  " + K_string[47]);
	        System.out.print("   TAT Y Tyr  " + K_string[61]);
	        System.out.print("   TGT C Cys  " + K_string[45] + "\n");
	        System.out.print("   TTC F Phe  " + K_string[55]);
	        System.out.print("   TCC S Ser  " + K_string[39]);
	        System.out.print("   TAC Y Tyr  " + K_string[53]);
	        System.out.print("   TGC C Cys  " + K_string[37] + "\n");
	        System.out.print("   TTA L Leu  " + K_string[62]);
	        System.out.print("   TCA S Ser  " + K_string[46]);
	        System.out.print("   TAA * Ter  " + K_string[60]);
	        System.out.print("   TGA * Ter  " + K_string[44] + "\n");
	        System.out.print("   TTG L Leu  " + K_string[54]);
	        System.out.print("   TCG S Ser  " + K_string[38]);
	        System.out.print("   TAG * Ter  " + K_string[52]);
	        System.out.print("   TGG W Trp  " + K_string[36] + "\n");
	
	        System.out.print("   CTT L Leu  " + K_string[31]);
	        System.out.print("   CCT P Pro  " + K_string[15]);
	        System.out.print("   CAT H His  " + K_string[29]);
	        System.out.print("   CGT R Arg  " + K_string[13] + "\n");
	        System.out.print("   CTC L Leu  " + K_string[23]);
	        System.out.print("   CCC P Pro  " + K_string[7]);
	        System.out.print("   CAC H His  " + K_string[21]);
	        System.out.print("   CGC R Arg  " + K_string[5] + "\n");
	        System.out.print("   CTA L Leu  " + K_string[30]);
	        System.out.print("   CCA P Pro  " + K_string[14]);
	        System.out.print("   CAA Q Gln  " + K_string[28]);
	        System.out.print("   CGA R Arg  " + K_string[12] + "\n");
	        System.out.print("   CTG L Leu  " + K_string[22]);
	        System.out.print("   CCG P Pro  " + K_string[6]);
	        System.out.print("   CAG Q Gln  " + K_string[20]);
	        System.out.print("   CGG R Arg  " + K_string[4] + "\n");
	
	        System.out.print("   ATT I Ile  " + K_string[59]);
	        System.out.print("   ACT T Thr  " + K_string[43]);
	        System.out.print("   AAT N Asn  " + K_string[57]);
	        System.out.print("   AGT S Ser  " + K_string[41] + "\n");
	        System.out.print("   ATC I Ile  " + K_string[51]);
	        System.out.print("   ACC T Thr  " + K_string[35]);
	        System.out.print("   AAC N Asn  " + K_string[49]);
	        System.out.print("   AGC S Ser  " + K_string[33] + "\n");
	        System.out.print("   ATA I Ile  " + K_string[58]);
	        System.out.print("   ACA T Thr  " + K_string[42]);
	        System.out.print("   AAA K Lys  " + K_string[56]);
	        System.out.print("   AGA R Arg  " + K_string[40] + "\n");
	        System.out.print("   ATG M Met  " + K_string[50]);
	        System.out.print("   ACG T Thr  " + K_string[34]);
	        System.out.print("   AAG K Lys  " + K_string[48]);
	        System.out.print("   AGG R Arg  " + K_string[32] + "\n");
	
	        System.out.print("   GTT V Val  " + K_string[27]);
	        System.out.print("   GCT A Ala  " + K_string[11]);
	        System.out.print("   GAT D Asp  " + K_string[25]);
	        System.out.print("   GGT G Gly  " + K_string[9] + "\n");
	        System.out.print("   GTC V Val  " + K_string[19]);
	        System.out.print("   GCC A Ala  " + K_string[3]);
	        System.out.print("   GAC D Asp  " + K_string[17]);
	        System.out.print("   GGC G Gly  " + K_string[1] + "\n");
	        System.out.print("   GTA V Val  " + K_string[26]);
	        System.out.print("   GCA A Ala  " + K_string[10]);
	        System.out.print("   GAA E Glu  " + K_string[24]);
	        System.out.print("   GGA G Gly  " + K_string[8] + "\n");
	        System.out.print("   GTG V Val  " + K_string[18]);
	        System.out.print("   GCG A Ala  " + K_string[2]);
	        System.out.print("   GAG E Glu  " + K_string[16]);
	        System.out.print("   GGG G Gly  " + K_string[0] + "\n");
		}else {
			String outfilename = "D:\\My Documents\\GThesis\\Res\\DNA_Beilin_out.txt";
			String strbuffer = new String();
			try {
				FileOutputStream outfile = new FileOutputStream(outfilename);
				DataOutputStream outdata = new DataOutputStream(outfile);
				
				strbuffer = new String("DNA_Beilin_out, start point of : " + start_point);
				outdata.writeUTF(strbuffer);
				
				byte k_string[] = new byte[K_dim];
				int k_loop = 0;
				int k = 0;
				int MAX = (int)Math.pow(4, K_dim);
				int loop_input = 0;
				long value = 0;
				
				outfile.write(13);
				
				while (k_loop < MAX) {
					outfile.write(10);
					for (k = 0; k < K_dim; k ++)
						k_string[k] = (byte)anti_coding(k_loop)[k];
					outfile.write(k_string);
					outfile.write(32);
					outfile.write(58);
					outfile.write(32);
					
					value = K_string[k_loop];
					strbuffer = (new Long(value)).toString();
					for (loop_input = 0; loop_input < strbuffer.length(); loop_input ++)
						outfile.write((byte)strbuffer.charAt(loop_input));
					
					outfile.write(13);
					k_loop ++;
				}
				outdata.close();
				outfile.close();
			}catch(IOException e) {
				throw e;
			}
		}
	}
	
    public static void main(String[] arg) throws IOException {
    	new DNA_Beilin_stat();
    }
}

