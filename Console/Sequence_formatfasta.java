package Console;

import java.io.*;

public class Sequence_formatfasta {

	/**
	 * @param args
	 */
	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		String infilename =
            "D:\\My Documents\\GThesis\\DatabaseBAK\\std_aa_1M_01.xml";
		String outfilename =
            "D:\\My Documents\\GThesis\\DatabaseBAK\\std_aa_1M_01_fasta.txt";
		
		try {
            FileInputStream infile = new FileInputStream(infilename);
            FileOutputStream outfile = new FileOutputStream(outfilename);
            
            DataOutputStream outdata = new DataOutputStream(outfile);
            String loop_str = new String();
            
            byte data[] = new byte[1];
            boolean start = false;
            long loop = 0;
            int i = 0;
            
            while (infile.read(data, 0, 1) != -1){
            	        	
            	if (data[0] == 61) {
            		loop ++;
            		outfile.write(62);
            		outfile.write(82);
            		outfile.write(65);
            		outfile.write(78);
            		outfile.write(68);
            		outfile.write(79);
            		outfile.write(77);
            		
            		loop_str = "_" + loop;
            		outdata.writeBytes(loop_str);
            		
	            	outfile.write(13);
	            	start = true;
            	}
            	if (start) {
            		outfile.write(10);
            		for (i = 0; i< 121; i++) {
            			infile.read(data, 0, 1);
            			if (data[0] > 64 && data[0] < 91)
            				outfile.write(data[0]);
            		}
            		outfile.write(13);
            		outfile.write(10);
            		start = false;
            	}
            }
            
            infile.close();
            outfile.close();
		}catch (IOException e) {
			throw e;
		}
		
		
	}

}
