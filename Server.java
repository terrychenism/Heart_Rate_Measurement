import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.ServerSocket;
import java.net.Socket;

public class Server {

	private int port = 9999;

	private ServerSocket ss;

	private Socket socket;

	/** Data input stream from Internet **/
	private DataInputStream dis;

	/** Data output stream to write to file **/
	private DataOutputStream dos;

	private String savePath = "D:\\";

	public void startServer() {

		// Percentage
		int passedlength = 0;

		// File length
		long length = 0;

		//Buffer size
		int bufferSize = 8192;

		// Buffer
		byte[] buf = new byte[bufferSize];

		try {

			ss = new ServerSocket(port);
			
			socket = ss.accept();
			
			System.out.println("-----connecting----");
			
			dis = new DataInputStream(new BufferedInputStream(
					socket.getInputStream()));

		} catch (IOException e) {

			e.printStackTrace();

		}

		try {
			// Read file name
			savePath += dis.readUTF();
			// Read file length
			length = dis.readLong();

		} catch (IOException e) {

			e.printStackTrace();

		}

		try {

			dos = new DataOutputStream(new BufferedOutputStream(
					new FileOutputStream(savePath)));

		} catch (FileNotFoundException e) {

			e.printStackTrace();

		}

		while (true) {

			int read = 0;

			if (dis != null) {

				try {

					read = dis.read(buf);

				} catch (IOException e) {

					e.printStackTrace();
				}
			}

			passedlength += read;

			if (read == -1) {

				break;

			}

			System.out.println("File Received:" + (passedlength * 100 / length) + "%\n");

			try {

				dos.write(buf, 0, read);

			} catch (IOException e) {

				e.printStackTrace();

			}

		}
		
		if (dos != null) {
			try {
			dos.close();
			} catch (IOException e) {
			e.printStackTrace();
		}
			dos = null;
		}
		
		if (dis != null) {
			try {
				dis.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			dis = null;
		}
		if (socket != null) {
			try {
				socket.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			socket = null;
		}
		
		if (ss != null) {
			try {
				ss.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
			ss = null;
		}
		System.out.println("Completed, file is located at:" + savePath + "\n");
		
	}
	public void start(String re) {
        try {
            ServerSocket serverSocket = new ServerSocket(port);
            while (true) {
            	socket = serverSocket.accept();
                try {
					System.out.println("-----connecting----");
					setMessage(re);
				} catch (Exception e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}finally{
		        	socket.close();
		        	serverSocket.close();
					System.out.println("close");
				}
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
	private void setMessage(String re) throws IOException{
    	String result=re.trim();
    	PrintWriter os=new PrintWriter(socket.getOutputStream());
        os.println(result);
        os.flush();
        System.out.println("-----Sending complete------");
        os.close();
        System.exit(0);
    }

    public static void readTxtFile(String filePath){
        try {
                String encoding="GBK";
                File file=new File(filePath);
                if(file.isFile() && file.exists()){ 
                    InputStreamReader read = new InputStreamReader(
                    new FileInputStream(file),encoding);
                    BufferedReader bufferedReader = new BufferedReader(read);
                    String lineTxt = null;
                    while((lineTxt = bufferedReader.readLine()) != null){
                        System.out.println(lineTxt);
                    }
                    read.close();
        }else{
            System.out.println("cannot received data");
        }
        } catch (Exception e) {
            System.out.println("error!!!");
            e.printStackTrace();
        }
     
    }

	public static void main(String[] args) {
		Server fileServer = new Server();
		String txtFilePath = "C:\\Users\\Terry\\Desktop\\exp.txt";

		String re=null;
		while(true){
		System.out.println("-----Ready to connect----");
		fileServer.startServer();
		System.out.println("-----Your heart rate----:");
		fileServer.readTxtFile(txtFilePath);
		BufferedReader in=new BufferedReader(new InputStreamReader(System.in)); 
		
		
		
		try {
			re=in.readLine().trim();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		System.out.println("-----Ready to connect----");
		fileServer.start(re);
		}
	}
}

