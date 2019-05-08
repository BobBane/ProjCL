/******************************************************************************
Copyright © 1999-2007, United States Government as represented by the Administrator for The National Aeronautics and Space Administration.  All Rights Reserved.
 *******************************************************************************/

package gov.nasa.gsfc.drl.opencl.projcl;

//import java.util.Arrays;

//import gov.nasa.gsfc.drl.h2gx.utilities.*;

/**
 * JNI interface class for ProjCL, a library that does proj.4 projections
 * using OpenCL
 *
 */
public class ProjCL {
    private static boolean initializedFlag = false;
    static {
	//System.err.println("Trying to initialize ProjCL...");
	try {
	    System.loadLibrary("projcl_jni_c");
	   
	    String initVal = null;
	    if ((initVal = init()) != null)
		throw new Exception(initVal);
	    
	    Runtime.getRuntime().addShutdownHook (new Thread()
		{
		    public void run() { quit(); }
		}
		);
	    initializedFlag = true;
		//System.err.println("ProjCL initialized");
	}
	catch (Throwable e) {
	    System.err.println("ProjCL initialization failed:");
	    System.err.println(e);
	    System.err.println("Continuing without ProjCL support.");
	    // e.printStackTrace();
	}
    }

    /**
     * True if the JNI part of the library is loaded and initialized
     */
    public static boolean initialized() {
	return initializedFlag;
    }

    /**
     * Initializes projcl package. Returns null if everything is OK,
     * String describing fault if not.
     */
    public static native String init();

    /**
     * Cleans up and releases openCL resources; called by shutdown
     * hook.
     */
    public static native void quit();

    /**
     * Looks up a proj.4 projection name and returns a non-zero integer
     * if ProjCL supports it, or zero if it does not.
     */
    public static native int
	getProjection(String p4name);

    /**
     * Performs forward ProjCL projection.  Returns NULL or a
     * string describing any problem.
     */
    public static native String
	projectForward (String p4name,
			double xyin[],
			double xyout[],
			double latOrigin, double lonOrigin,
			double standardParallel1,
			double standardParallel2
			);

    /**
     * Performs forward ProjCL projection.  Returns NULL or a
     * string describing any problem.
     */
    public static native String
	projectForward2 (String p4name,
			double xin[], double yin[],
			double xout[], double yout[],
			double latOrigin, double lonOrigin,
			double standardParallel1,
			double standardParallel2
			);



    /**
     * Trivial test
     */
    public static void main(String[] args) {
	String initVal = null;
	if (initialized()) {
	    System.err.println("ProjCL initialized");
	}
	else
	    System.err.println("ProjCL NOT initialized - " + initVal);
	System.err.println("stere: " + getProjection("stere"));
	System.err.println("sterea: " + getProjection("sterea"));
    }
}
