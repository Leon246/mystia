package Visualize;

import java.io.*;
import java.awt.image.*;
import javax.imageio.*;
import java.awt.*;

public class Theoretical_fractal {
	public int K_dim = 10;
	
	public int seq_exists(int X, int Y, String seq) {
		char k_string[] = new char[K_dim];
		char k_temp = 'U';
		int line = (int)Math.pow(2, K_dim);
		
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
		String str = new String(k_string);
		if (str.contains(seq)) {
			return 1;
		}
		else {
			return 0;
		}
	}
	
	public void Theoretical_fractal_single(String missing_seq) throws IOException{
		String spe = missing_seq;
		if (spe.length() < K_dim) {

	        int i = 0, j = 0;
			int x = 1, y = 1;
			int MAX = (int)Math.pow(2, K_dim);
			int a = 1024 / MAX;
			BufferedImage bufferedImage = new BufferedImage(1026, 1026, BufferedImage.TYPE_INT_RGB);
			    
		    // Create a graphics contents on the buffered image
		    Graphics2D g2d = bufferedImage.createGraphics();
		    
			for (i = 0; i < MAX; i ++) {
				for (j = 0; j < MAX; j ++) {
					if (seq_exists(i, j, spe) == 1) {
						g2d.setColor(Color.black);
					}else {
						g2d.setColor(Color.gray);
					}
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
	            File file = new File("D:\\My Documents\\GThesis\\Res\\image\\" + spe + "_" + K_dim + "_lack.png");
	            ImageIO.write(bufferedImage, "png", file);
	            System.out.print("Fractal Printed\n");
	        } catch (IOException e) {
	        	throw e;
	        }
		}else {
			System.out.print("special sequence cannot be too long!\n");
			System.exit(0);
		}	
	}
	
	public void Theoretical_fractal_muti() throws IOException{
		String spe = new String("TA");
		String spe_2 = new String("TATA");
		String spe_3 = new String("CTAA");
		String spe_4 = new String("CTAG");
		if (spe.length() < K_dim) {
			int width = 1026;
			int height = 1026;

	        int i = 0, j = 0;
			int x = 1, y = 1;
			int MAX = (int)Math.pow(2, K_dim);
			int a = 1024 / MAX;
			int color = 0;
			BufferedImage bufferedImage = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
			    
		    // Create a graphics contents on the buffered image
		    Graphics2D g2d = bufferedImage.createGraphics();
		    int RGB = 0;
		    
			for (i = 0; i < MAX; i ++) {
				for (j = 0; j < MAX; j ++) {
					color = seq_exists(i, j, spe) + seq_exists(i, j, spe_2) + seq_exists(i, j, spe_3) + seq_exists(i, j, spe_4);
					if (color == 0) {
						RGB = 255;
					}else if (color == 1) {
						RGB = 64;
					}else if (color == 2) {
						RGB = 32;
					}else if (color == 3) {
						RGB = 16;
					}else if (color == 4) {
						RGB = 0;
					}
					
					
					x = i * a + 1;
					y = j * a + 1;
					Color co = new Color(RGB, RGB, RGB);
					g2d.setColor(co);
					g2d.fillRect(x, y, x + a, y + a);
				}
			}
	        // Graphics context no longer needed so dispose it
	        g2d.dispose();
	        
	        // Write generated image to a file
	        try {
	            // Save as PNG
	            File file = new File("D:\\My Documents\\GThesis\\Res\\image\\muti_lack.png");
	            ImageIO.write(bufferedImage, "png", file);
	            System.out.print("Fractal Printed\n");
	        } catch (IOException e) {
	        	throw e;
	        }
		}else {
			System.out.print("special sequence cannot be too long!\n");
			System.exit(0);
		}	
	}
	
	public Theoretical_fractal() throws IOException{
		//Theoretical_fractal_single("TA");
		Theoretical_fractal_muti();
	}
	
	public static void main(String arg[]) throws IOException{

		/*
		char na[] = new char[4];
		na[0] = 'G';
		na[1] = 'C';
		na[2] = 'A';
		na[3] = 'T';
		
		int na_serial = 0;
		int i = 0;
		int seq_length = 1;
		while (seq_length < 5) {
			for (na_serial = 0; na_serial < 4; na_serial ++) {
				String lack_seq = new String();
				for (i = 0; i < seq_length; i++) {
					lack_seq += na[na_serial];
				}
				System.out.print(lack_seq + "\n");
				new Theoretical_fractal(lack_seq);
			}
			seq_length ++;
		}*/
		new Theoretical_fractal();
	}
}