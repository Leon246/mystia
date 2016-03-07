package Visualize;

import java.io.*;
import java.awt.image.*;
import javax.imageio.*;
import java.awt.*;

public class Pic_Generator {
	
//	set K_string dimension
	public int K_dim = 10;
	public int start_point = 0;
	
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
	
	public void DNA_Beilin_stat(String filename) throws IOException {

        //Define K_string, need not to be reset to zero
        int K = K_dim;
        int k = 0;
        int n = 0;
        K_string = new long[(int) Math.pow(4, (double)K)];

        try {
            FileInputStream infile = new FileInputStream(filename);
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

            //close file
            indata.close();
            infile.close();
        } catch (IOException e) {
            throw e;
        }
	}
	
	public void DNA_normal_stat(String filename) throws IOException {

	    //Define K_string, need not to be reset to zero
	    //K defines the K_string matrix, significant to the most
	    int K = K_dim;
	    int k = 0;
	    int n = 0;
	    
	    String strbuffer = new String();
	
	    K_string = new long[(int) Math.pow(4, (double) K)];
	    
	    try {
	        FileInputStream infile = new FileInputStream(filename);
	        DataInputStream indata = new DataInputStream(infile);
	
	        int X = 0, Y = 0, value = 0;
	        int half = (int) Math.pow(2., (double) (K - 1)), total = half * 2;
	        byte data[] = new byte[1];
	        byte store[] = new byte[K_dim];
	        int store_loop = start_point;
	        if (store_loop >= K_dim) {
	        	System.out.print("Start point must be no more than K_dim! \n");
	        	System.exit(2);
	        }
	        char seq[] = new char[K];
	        
	        boolean firstline = true;
	        int nbyteread = indata.read(data, 0, 1);
	
	        while (nbyteread != -1) {
	        	
	        	//be very careful with the .fna format
	        	//jump over the first line only
	        	if (firstline)
	        		nbyteread = indata.read(data, 0, 1);
	        	
	        	if (data[0] == 10 && firstline) {
	        		firstline = false;
	    	        System.out.print("First line passed! \n");
		        }else if (!firstline){
		            if ((data[0] == 10 || data[0] == 13) && !firstline) {
		            	nbyteread = indata.read(data, 0, 1);
		                continue;
		            } else {
		                if (k == K) {	                	
		                	//read now
		                    while (k > 0) {
		                        //X-axis of sequence
		                        value = X_coding(seq[k - 1]);
		                        X += (int) Math.pow(2, K - k) * value;
		                        //Y-axis of sequence
		                        value = Y_coding(seq[k - 1]);
		                        Y += (int) Math.pow(2, K - k) * value;
		                        k--;
		                    }
		                    if (k == 0) {
		                        K_string[X + Y * total]++;
		                        X = 0;
		                        Y = 0;
		                        n += K_dim;
		                        //temp:
		                        //strbuffer = new String(seq);
		                        //System.out.print(strbuffer + "\n");
		                    } else {
		                        System.out.print("k != 0, Error");
		                    }
		                } else if (k < K) {
		                	//set start point
		                	while (store_loop > 0) {
		                		store[start_point - store_loop] = data[0];
		                		indata.read(data, 0, 1);
		                		store_loop --;
		                	}
		                    seq[k] = (char) data[0];
		                    k ++;
		                    if (!firstline)
		                    	nbyteread = indata.read(data, 0, 1);
		                } else {
		                    System.out.print("Unexceptable Error");
		                }
		            }
	        	}else {
	        		continue;
	        	}
	        }
	        //read as loop
	       	store_loop = 0;
	       	int store_remains = start_point;
	       	while (k < K_dim && store_remains > 0) {
	       		seq[k] = (char) store[store_loop];
	       		store_loop ++;
	       		k ++;
	       		store_remains --;
	       	}
	        if (k == K) {	                	
            	//read now
                while (k > 0) {
                    //X-axis of sequence
                    value = X_coding(seq[k - 1]);
                    X += (int) Math.pow(2, K - k) * value;
                    //Y-axis of sequence
                    value = Y_coding(seq[k - 1]);
                    Y += (int) Math.pow(2, K - k) * value;
                    k--;
                }
                if (k == 0) {
                    K_string[X + Y * total]++;
                    X = 0;
                    Y = 0;
                    n += K_dim;
                    //temp:
                    strbuffer = new String(seq);
                    System.out.print(strbuffer + "\n");
                } else {
                    System.out.print("k != 0, Error");
                }
            }
			
            //close file
            indata.close();
            infile.close();
	    }catch(IOException e) {
	    	throw e;
	    }
	}
	
	public Pic_Generator(String filename, String imagename) throws IOException{
		String image_name = new String();

		DNA_Beilin_stat(filename);
//		DNA_normal_stat(filename);
		
        int i = 0, j = 0;
		int x = 1, y = 1;
		int R = 0, G = 0, B = 0;
		int K_value = 0;
		int MAX = (int)Math.pow(2, K_dim);
		int a = 1024 / MAX;
		
		int color_set = 20;
		
//		for (color_set = 3; color_set < 42; color_set += 2) {

			BufferedImage bufferedImage = new BufferedImage(1026, 1026, BufferedImage.TYPE_INT_RGB);
		    
	        // Create a graphics contents on the buffered image
	        Graphics2D g2d = bufferedImage.createGraphics();

			for (i = 0; i < MAX; i ++) {
				for (j = 0; j < MAX; j ++) {
					K_value = (int)(K_string[i + j * MAX]);
					K_value = K_value * color_set;
					if (K_value < 256) {
						R = K_value / 3;
						G = K_value / 1;
						B = K_value / 3;
					}else if (K_value >= 256){
						R = 255;
						G = 255;
						B = 255;
					}else {
						System.out.print("Error in setting color");
						System.exit(0);
					}
					
					Color co = new Color(R, G, B);
					g2d.setColor(co);
					x = i * a + 1;
					y = j * a + 1;
					g2d.fillRect(x, y, x + a, y + a);
				}
			}
	        // Graphics context no longer needed so dispose it
	        g2d.dispose();
	        
	        // Write generated image to a file
	        try {
	            // Save as PNG
	        	image_name = imagename + "_" +K_dim + "_" + start_point + "_" +color_set + ".png";
	            File file = new File(image_name);
	            ImageIO.write(bufferedImage, "png", file);
	        } catch (IOException e) {
	        	throw e;
	        }
//		}
	}
	
	public static void main(String arg[]) throws IOException{
		int file_num = 0;
		String temp = new String();
	    String search = 
	    	"D:\\My Documents\\GThesis\\fna_data\\NC_";
	    String result =
	    	"D:\\res\\NC_";
	    String filename = new String();
		String imagename = new String();
		for (file_num = 913; file_num < 914; file_num ++) {
			temp = "0" + file_num;
			while (temp.length() < 6) {
				temp = "0" + temp;
			}
			filename = search + temp + ".fna";
			File f = new File(filename);
			if (f.exists()) {
				imagename = result + temp;
				new Pic_Generator(filename, imagename);
			}
		}
	}
}