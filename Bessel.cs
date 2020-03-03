using System;
using System.Numerics;

namespace ComplexBessel
{
	/// <summary>
	/// Contains integer order Bessel's Functions of real and complex argument.
	/// \todo J_n, Y_n for complex z, n = 0,...; I_n, K_n for complex argument, n = 1, 2, ... .
	/// </summary>
	public static class Bessel
	{
		#region Constants
		// γ (Euler-Mascheroni)
		public const double Gamma			= 0.577215664901532860606512;
		public const double OneOverSqrt2 	= 0.7071067811865475244;
		public const double FineAccuracy	= 1e-9;
		public const double MachineAccuracy	= 1e-14;
		#endregion

		#region Functions J_n(x)
		/// <summary>
		/// Definition: Abramowitz and Steagun 9.4.1 (http://people.math.sfu.ca/~cbm/aands/page_358.htm)
		/// With |x| &lt;= 8 Algorithm JZERO 5842 from J.F. Hart, Computer Approximations is used.
		/// With |x| &gt; 8 Asymptotic approximation with coefficients borrowed from Cephes Library.
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value of J0(x) with accuracy ~ 1e-11.</returns>
		public static double J0(double x)
		{
			if ((Math.Abs(x)) <= 8.0)
			{
				return J0Hart5342(x);
			}
			else
			{
				x = Math.Abs(x);

				double pzero = 0;
				double qzero = 0;

				J0Asymptotic(x, ref pzero, ref qzero);
				double nn = x - Math.PI / 4;
				double dResult = Math.Sqrt(2 / Math.PI / x) * (pzero * Math.Cos(nn) - qzero * Math.Sin(nn));

				return dResult;
			}
		}

		/// <summary>
		/// Definition: Abramowitz and Steagun 9.4.1 (http://people.math.sfu.ca/~cbm/aands/page_358.htm)
		/// With |x| &lt;= 8 Algorithm JZERO 6042 from J.F. Hart, Computer Approximations is used.
		/// With |x| &gt; 8 Asymptotic approximation with coefficients borrowed from Cephes Library.
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value of J1(x) with accuracy ~ 1e-11.</returns>
		/// <remarks>
		///		Since values of J1 are given with accuracy 10^-8 in Abramowitz and Steagun, this function is only tested to this accuracy.
		///		Assertion that it delivers accuracy of 10^-11 is based on Hart's tests.
		///	</remarks>
		public static double J1(double x)
		{
			if ((Math.Abs(x)) <= 8.0)
			{
				return x * J1Hart6042(x);
			}
			else
			{
				double x1 = Math.Abs(x);

				double pzero = 0;
				double qzero = 0;

				J1Asymptotic(x1, ref pzero, ref qzero);

				double nn = x1 - 3 * Math.PI / 4;

				double dResult = Math.Sqrt(2 / Math.PI / x1) * (pzero * Math.Cos(nn) - qzero * Math.Sin(nn));

				if (x < 0)
				{
					dResult = -dResult;
				}

				return dResult;
			}
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="n"></param>
		/// <param name="x"></param>
		/// <returns></returns>
		public static double Jn(int n, double x)
		{
			switch (n)
			{
				case 0:
					return J0(x);

				case 1:
					return J1(x);

				default:
					double j_n_1	= J0(x);
					double j_n		= J1(x);

					for (int i = 1; i < n; i++)
					{
						double temp	= j_n;
						j_n			= (2 * i / x) * j_n - j_n_1;
						j_n_1		= temp;
					}

					return j_n;
			}
		}
		#endregion

		#region Functions Y_n(x)
		/// <summary>
		/// 
		/// </summary>
		/// <param name="x"></param>
		/// <returns></returns>
		public static double Y0(double x)
		{
			if ((Math.Abs(x)) <= 8.0)
			{
				return Y0Hart6236(x) + (2.0 / Math.PI) * J0(x) * Math.Log(x);
			}
			else
			{
				double pzero = 0;
				double qzero = 0;

				Y0Asymptotic(x, ref pzero, ref qzero);

				double nn = x - Math.PI / 4;
				return Math.Sqrt(2 / Math.PI / x) * (pzero * Math.Sin(nn) + qzero * Math.Cos(nn));
			}
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="x"></param>
		/// <returns></returns>
		public static double Y1(double x)
		{
			if ((Math.Abs(x)) <= 8.0)
			{
				return Y1Hart6438(x) + (2.0 / Math.PI) * (J1(x) * Math.Log(x) - 1.0 / x);
			}
			else
			{
				x				= Math.Abs(x);
				double pzero	= 0;
				double qzero	= 0;

				J1Asymptotic(x, ref pzero, ref qzero);
				double nn		= x - 3 * Math.PI / 4;
				return Math.Sqrt(2 / Math.PI / x) * (pzero * Math.Sin(nn) + qzero * Math.Cos(nn));
			}
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="n"></param>
		/// <param name="x"></param>
		/// <returns></returns>
		public static double Yn(int n, double x)
		{
			switch (n)
			{
				case 0:
					return Y0(x);

				case 1:
					return Y1(x);

				default:
					double j_n_1	= Y0(x);
					double j_n		= Y1(x);

					for (int i = 1; i < n; i++)
					{
						double temp	= j_n;
						j_n			= (2 * i / x) * j_n - j_n_1;
						j_n_1		= temp;
					}

					return j_n;
			}
		}
		#endregion

		#region Functions I_n(x)
		/// <summary>
		/// Modified Bessel function of first kind, order 0 for real argument.
		/// Definition: Abramowitz and Steagun 9.6          (http://www.math.sfu.ca/%7Ecbm/aands/page_374.htm)
		/// Source:     Abramowitz and Steagun 9.8.1, 9.8.2 (http://www.math.sfu.ca/%7Ecbm/aands/page_378.htm)
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value</returns>
		public static double I0(double x)
		{
			const double a375 = 3.75;

			double[] a = new double[]
			{
				1.0, 		3.5156229, 3.0899424,
				1.2067492,	0.2659732, 0.0360768, 0.0045813
			};

			double t = x / a375;
			double s;

			if (Math.Abs(t) <= 1.0)
			{
				t = t * t;
				s = a[6] * t + a[5];
				s = s * t + a[4];
				s = s * t + a[3];
				s = s * t + a[2];
				s = s * t + a[1];
				s = s * t + a[0];
			}
			else
			{
				t = 1.0 / t;

				double[] b = new double[]
				{
					0.39894228,  0.01328592, 0.00225319, -0.00157565,
					0.00916281, -0.02057706, 0.02635537, -0.01647633,
					0.00392377
				};

				s = b[8] * t + b[7];
				s = s * t + b[6];
				s = s * t + b[5];
				s = s * t + b[4];
				s = s * t + b[3];
				s = s * t + b[2];
				s = s * t + b[1];
				s = s * t + b[0];
				s = s * Math.Exp(x) / Math.Sqrt(x);
			}

			return s;
		}

		/// <summary>
		/// Modified Bessel function of first kind, order 1 for real argument.
		/// Definition: Abramowitz and Steagun 9.6				(http://www.math.sfu.ca/%7Ecbm/aands/page_374.htm)
		/// Source:     Abramowitz and Steagun 9.8.3, 9.8.4		(http://www.math.sfu.ca/%7Ecbm/aands/page_378.htm)
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value</returns>
		public static double I1(double x)
		{
			const double a375 = 3.75;

			double[] a = new double[]
			{
				0.5, 		 0.87890594, 0.51498869, 0.15084934,
				0.026858733, 0.00301532, 0.00032411
			};

			double t = x / a375;
			double s;

			if (Math.Abs(t) <= 1.0)
			{
				t = t * t;
				s = a[6] * t + a[5];
				s = s * t + a[4];
				s = s * t + a[3];
				s = s * t + a[2];
				s = s * t + a[1];
				s = s * t + a[0];
				s = s * x;
			}
			else
			{
				t = 1.0 / t;

				double[] b = new double[]
				{
					0.39894228, -0.03988024, -0.00362018,
					0.00163801, -0.01031555,  0.02282967,
					-0.02895312, 0.01787654, -0.00420059
				};

				s = b[8] * t + b[7];
				s = s * t + b[6];
				s = s * t + b[5];
				s = s * t + b[4];
				s = s * t + b[3];
				s = s * t + b[2];
				s = s * t + b[1];
				s = s * t + b[0];
				s = s * Math.Exp(x) / Math.Sqrt(x);
			}

			return s;
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="n"></param>
		/// <param name="x"></param>
		/// <returns></returns>
		public static double In(int n, double x)
		{
			switch (n)
			{
				case 0:
					return I0(x);

				case 1:
					return I1(x);

				default:
					double j_n_1	= I0(x);
					double j_n		= I1(x);

					for (int i = 1; i < n; i++)
					{
						double temp	= j_n;
						j_n			= -(2 * i / x) * j_n + j_n_1;
						j_n_1		= temp;
					}

					return j_n;
			}
		}
		#endregion

		#region Functions K_n(x)
		/// <summary>
		/// Modified Bessel function of second kind, order 0 (Macdonald function of order 0) for real argument.
		/// Definition: Abramowitz and Steagun 9.6          (http://www.math.sfu.ca/%7Ecbm/aands/page_374.htm)
		/// Source:     Abramowitz and Steagun 9.8.5, 9.8.6 (http://www.math.sfu.ca/%7Ecbm/aands/page_379.htm)
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value</returns>
		public static double K0(double x)
		{
			double t;
			double s;
			if ((0.0 < x) && (x <= 2.0))
			{
				double[] a = new double[]
				{
					-0.57721566, 0.42278420, 0.23069756, 0.03488590,
					0.00262698,  0.00010750, 0.00000740
				};

				t = 0.5 * x;
				t = t * t;
				s = a[6] * t + a[5];
				s = s * t + a[4];
				s = s * t + a[3];
				s = s * t + a[2];
				s = s * t + a[1];
				s = s * t + a[0];

				return s - I0(x) * Math.Log(0.5 * x);
			}
			else
			{
				double[] b = new double[]
				{
					1.25331414, -0.07832358, 0.02189568, -0.01062446,
					0.00587872, -0.00251540, 0.00053208
				};

				t = 2.0 / x;
				s = b[6] * t + b[5];
				s = s * t + b[4];
				s = s * t + b[3];
				s = s * t + b[2];
				s = s * t + b[1];
				s = s * t + b[0];

				return s * Math.Exp(-x) / Math.Sqrt(x);
			}
		}

		/// <summary>
		/// Modified Bessel function of second kind, order 1 (Macdonald function of order 1) for real argument.
		/// Definition: Abramowitz and Steagun 9.6				(http://www.math.sfu.ca/%7Ecbm/aands/page_374.htm)
		/// Source:     Abramowitz and Steagun 9.8.7, 9.8.8		(http://www.math.sfu.ca/%7Ecbm/aands/page_379.htm)
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value</returns>
		public static double K1(double x)
		{
			double t;
			double s;
			if ((0.0 < x) && (x <= 2.0))
			{
				double[] a = new double[]
				{
					1.0, 		  0.15443144, -0.67278579, -0.18156897,
					-0.01919402, -0.00110404, -0.00004686
				};

				t = 0.5 * x;
				t = t * t;
				s = a[6] * t + a[5];
				s = s * t + a[4];
				s = s * t + a[3];
				s = s * t + a[2];
				s = s * t + a[1];
				s = s * t + a[0];

				return I1(x) * Math.Log(0.5 * x) + s / x;
			}

			else
			{
				double[] b = new double[]
				{
					1.25331414, 0.23498619, -0.03655620, 0.01504268,
					-0.00780353, 0.00325614, -0.00068245
				};

				t = 2.0 / x;
				s = b[6] * t + b[5];
				s = s * t + b[4];
				s = s * t + b[3];
				s = s * t + b[2];
				s = s * t + b[1];
				s = s * t + b[0];

				return s * Math.Exp(-x) / Math.Sqrt(x);
			}
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="n"></param>
		/// <param name="x"></param>
		/// <returns></returns>
		public static double Kn(int n, double x)
		{
			switch (n)
			{
				case 0:
					return K0(x);

				case 1:
					return K1(x);

				default:
					double j_n_1	= K0(x);
					double j_n		= K1(x);

					for (int i = 1; i < n; i++)
					{
						double temp	= j_n;
						j_n			= (2 * i / x) * j_n + j_n_1;
						j_n_1		= temp;
					}

					return j_n;
			}
		}
		#endregion

		#region Kelvin's Functions
		/// <summary>
		/// Real Kelvin function.
		/// Definition: Abramowitz and Steagun 9.9				(http://www.math.sfu.ca/~cbm/aands/page_379.htm)
		/// Source:     Abramowitz and Steagun 9.11.1, 9.11.9	(http://www.math.sfu.ca/%7Ecbm/aands/page_384.htm)
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value</returns>
		public static double Ber(double x)
		{
			if (Math.Abs(x) <= 8.0)
			{
				double t = x / 8;

				t *= t;
				t *= t;	// ^4

				double[] b = new double[]
				{
					1.0,		-64.0,		113.77777774,	-32.36345652,
					2.64191397, -0.8349609, 0.00122552,		-0.00000901
				};

				double s;
				s = b[7] * t + b[6];
				s = s * t + b[5];
				s = s * t + b[4];
				s = s * t + b[3];
				s = s * t + b[2];
				s = s * t + b[1];
				s = s * t + b[0];

				return s;
			}
			else if (x >= 8.0)
			{
				return (1.0 / (2.0 * Math.PI * x)) *
					   Math.Exp(OneOverSqrt2 * x + ThetaReal(x)) *
					   Math.Cos(OneOverSqrt2 * x + ThetaImagine(x)) -
					   (1.0 / Math.PI) * Kei(x);
			}
			else
			{
				return -Double.MaxValue;
			}
		}

		/// <summary>
		/// Imaginary Kelvin function.
		/// Definition: Abramowitz and Steagun 9.9				(http://www.math.sfu.ca/~cbm/aands/page_379.htm)
		/// Source:     Abramowitz and Steagun 9.11.1, 9.11.9	(http://www.math.sfu.ca/%7Ecbm/aands/page_384.htm)
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value</returns>
		public static double Bei(double x)
		{
			if (Math.Abs(x) <= 8.0)
			{
				double t = x / 8;

				t *= t;
				t *= t;	// ^4

				double[] b = new double[]
				{
					16.0,			-113.77777774,	72.81777742,
					-10.56765779,	0.52185615,		-0.01103667,	0.00011346
				};

				double s;
				s = b[6] * t + b[5];
				s = s * t + b[4];
				s = s * t + b[3];
				s = s * t + b[2];
				s = s * t + b[1];
				s = s * t + b[0];

				return s * x * x / 64;
			}
			else if (x >= 8.0)
			{
				return (1.0 / (2.0 * Math.PI * x)) *
					   Math.Exp(OneOverSqrt2 * x + ThetaReal(x)) *
					   Math.Sin(OneOverSqrt2 * x + ThetaImagine(x)) +
					   (1.0 / Math.PI) * Ker(x);
			}
			else
			{
				return -Double.MaxValue;
			}
		}

		/// <summary>
		/// Real Kelvin function.
		/// Definition: Abramowitz and Steagun 9.9				(http://www.math.sfu.ca/~cbm/aands/page_379.htm)
		/// Source:     Abramowitz and Steagun 9.11.1, 9.11.9	(http://www.math.sfu.ca/%7Ecbm/aands/page_384.htm)
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value</returns>
		public static double Ker(double x)
		{
			if (0.0 < x && x <= 8.0)
			{
				double t = x / 8;

				t *= t;
				t *= t;	// ^4

				double[] b = new double[]
				{
					-0.57721566,	-59.05819744,	171.36272133,	-60.60977451,
					5.65539121,		-0.19636347,	0.0309699,		-0.00002458
				};

				double s;
				s = b[7] * t + b[6];
				s = s * t + b[5];
				s = s * t + b[4];
				s = s * t + b[3];
				s = s * t + b[2];
				s = s * t + b[1];
				s = s * t + b[0];

				return s + (Math.PI / 4) * Bei(x) - Ber(x) * Math.Log(x / 2);
			}
			else if (x >= 8.0)
			{
				return Math.Sqrt(Math.PI / (2.0 * x)) *
					   Math.Exp(ThetaReal(-x) - x * OneOverSqrt2) *
					   Math.Cos(ThetaImagine(-x) - x * OneOverSqrt2);
			}
			else
			{
				return -Double.MaxValue;
			}
		}

		/// <summary>
		/// Imaginary Kelvin function.
		/// Definition: Abramowitz and Steagun 9.9				(http://www.math.sfu.ca/~cbm/aands/page_379.htm)
		/// Source:     Abramowitz and Steagun 9.11.1, 9.11.9	(http://www.math.sfu.ca/%7Ecbm/aands/page_384.htm)
		/// </summary>
		/// <param name="x">Argument</param>
		/// <returns>Value</returns>
		public static double Kei(double x)
		{
			if (0.0 < x && x <= 8.0)
			{
				double t = x / 8;

				t *= t;
				t *= t;	// ^4

				double[] b = new double[]
				{
					6.76454936,		-142.91827687,	124.23569650,
					-21.30060904,	1.17509064,		-0.02695875,	0.00029532
				};

				double s;
				s = b[6] * t + b[5];
				s = s * t + b[4];
				s = s * t + b[3];
				s = s * t + b[2];
				s = s * t + b[1];
				s = s * t + b[0];

				return s * x * x / 64 - (Math.PI / 4) * Ber(x) - Bei(x) * Math.Log(x / 2);
			}
			else if (x >= 8.0)
			{
				return Math.Sqrt(Math.PI / (2.0 * x)) *
					   Math.Exp(ThetaReal(-x) - x * OneOverSqrt2) *
					   Math.Sin(ThetaImagine(-x) - x * OneOverSqrt2);
			}
			else
			{
				return -Double.MaxValue;
			}
		}
		#endregion

		#region Complex Bessel		
		/// <summary>
		/// Computes I0(z) for a complex argument z.
		/// </summary>
		/// <param name="z">The complex argument.</param>
		/// <returns>The value of I0(z).</returns>
		public static Complex I0(Complex z)
		{
			if (z.Magnitude < _maxSeriesArgumentI0)
			{
				return ComplexI0BySeries(z);
			}
			else
			{
				return ComplexI0ByAsymptotic(z);
			}
		}

		/// <summary>
		/// Computes K0(z) for a complex argument z.
		/// </summary>
		/// <param name="z">The complex argument.</param>
		/// <returns>The value of K0(z).</returns>
		public static Complex K0(Complex z)
		{
			if (z.Magnitude < FineAccuracy)
			{
				throw new ArgumentException($"The argument of K0() is too small: {z}");
			}
			
			if (z.Magnitude < _maxSeriesArgumentK0)
			{
				return ComplexK0BySeries(z);
			}
			else
			{
				return ComplexK0ByAsymptotic(z);
			}
		}
		#endregion

		#region Private Auxiliary
		/// <summary>
		/// Implements calculation of J0(x) for |x| &lt;= 8 by means of rational polynomial approximation.
		/// (J.F. Hart, Computer Approximations, 1975. ISBN-13 978-0882756424, Algorithm 5842).
		/// </summary>
		/// <param name="x">Argument (|x| &lt; 8).</param>
		/// <returns>J0(x) with accuracy of approximately 1e-11.</returns>
		private static double J0Hart5342(double x)
		{
			double[] P	=	new double[]
										{
											+0.9005822861680792429196861076E14,
											-0.211835035258657051110887761E14,
											+0.108444107395388390805876161E13,
											-0.2075034062182266636524500445E11,
											+0.1774486522551782526401959939E9,
											-0.6921661031283507141299515818E6,
											+0.1012719509763049027304358973E4
										};

			double[] Q = new double[]
										{		
											+0.9005822861693449907844140699E14,
											+0.1331053627607928163120115803E13,
											+0.1004465946632439368187662601E11,
											+0.5088386102539175620522587556E8,
											+0.1896688971945018803547991105E6,
											+0.5244554464059932278280464543E3,
											1.0
										};

			return PolynomialValue(P, x * x) / PolynomialValue(Q, x * x);
		}
	
		/// <summary>
		/// Asymptotic Hankel decomposition with coefficients borrowed from Cephes library.
		/// </summary>
		/// <param name="x"></param>
		/// <param name="pzero"></param>
		/// <param name="qzero"></param>
		private static void J0Asymptotic(double x, ref double pzero, ref double qzero)
		{
			double z = 64.0 / (x * x);

			double p2 = 0.0;
			p2 = 2485.271928957404011288128951 + z * p2;
			p2 = 153982.6532623911470917825993 + z * p2;
			p2 = 2016135.283049983642487182349 + z * p2;
			p2 = 8413041.456550439208464315611 + z * p2;
			p2 = 12332384.76817638145232406055 + z * p2;
			p2 = 5393485.083869438325262122897 + z * p2;

			double q2 = 1.0;
			q2 = 2615.700736920839685159081813 + z * q2;
			q2 = 156001.7276940030940592769933 + z * q2;
			q2 = 2025066.801570134013891035236 + z * q2;
			q2 = 8426449.050629797331554404810 + z * q2;
			q2 = 12338310.22786324960844856182 + z * q2;
			q2 = 5393485.083869438325560444960 + z * q2;

			double p3 = -0.0;
			p3 = -4.887199395841261531199129300 + z * p3;
			p3 = -226.2630641933704113967255053 + z * p3;
			p3 = -2365.956170779108192723612816 + z * p3;
			p3 = -8239.066313485606568803548860 + z * p3;
			p3 = -10381.41698748464093880530341 + z * p3;
			p3 = -3984.617357595222463506790588 + z * p3;

			double q3 = 1.0;
			q3 = 408.7714673983499223402830260 + z * q3;
			q3 = 15704.89191515395519392882766 + z * q3;
			q3 = 156021.3206679291652539287109 + z * q3;
			q3 = 533291.3634216897168722255057 + z * q3;
			q3 = 666745.4239319826986004038103 + z * q3;
			q3 = 255015.5108860942382983170882 + z * q3;

			pzero = p2 / q2;
			qzero = 8 * p3 / q3 / x;
		}

		
		/// <summary>
		/// Implements calculation of J1(x) for |x| &lt;= 8 by means of rational polynomial approximation.
		/// (J.F. Hart, Computer Approximations, 1975. ISBN-13 978-0882756424, Algorithm 5842).
		/// </summary>
		/// <param name="x">Argument (|x| &lt; 8).</param>
		/// <returns>J1(x) with accuracy of approximately 1e-11.</returns>
		private static double J1Hart6042(double x)
		{
			double[] P = new double[]
										{
											+0.20401207733056935676404989E8,
											-0.24090488178461116078620946E7,					
											+0.89022956297134646818640971E5,					
											-0.152882002373665370815863364E4,					
											+0.14404277856294376706342034E2,					
											-0.804888363008675152104385216E-1,					
											+0.27202107121267732717690412E-3,					
											-0.5278642496908269100812841E-6,						
											+0.465970885798957648299543E-9				
										};

			double[] Q = new double[]
										{		
											+0.40802415466094926655209805E8,					
											+0.28220429766095254592028E6,						
											+0.808869176790579484514931E3,						
											+1.0
										};

			return PolynomialValue(P, x * x) / PolynomialValue(Q, x * x);
		}

		/// <summary>
		/// Asymptotic Hankel decomposition with coefficients borrowed from Cephes library.
		/// </summary>
		/// <param name="x"></param>
		/// <param name="pzero"></param>
		/// <param name="qzero"></param>
		private static void J1Asymptotic(double x, ref double pzero, ref double qzero)
		{
			double xsq = 64.0 / (x * x);
			double p2 = -1611.616644324610116477412898;
			p2 = -109824.0554345934672737413139 + xsq * p2;
			p2 = -1523529.351181137383255105722 + xsq * p2;
			p2 = -6603373.248364939109255245434 + xsq * p2;
			p2 = -9942246.505077641195658377899 + xsq * p2;
			p2 = -4435757.816794127857114720794 + xsq * p2;

			double q2 = 1.0;
			q2 = -1455.009440190496182453565068 + xsq * q2;
			q2 = -107263.8599110382011903063867 + xsq * q2;
			q2 = -1511809.506634160881644546358 + xsq * q2;
			q2 = -6585339.479723087072826915069 + xsq * q2;
			q2 = -9934124.389934585658967556309 + xsq * q2;
			q2 = -4435757.816794127856828016962 + xsq * q2;

			double p3 = 35.26513384663603218592175580;
			p3 = 1706.375429020768002061283546 + xsq * p3;
			p3 = 18494.26287322386679652009819 + xsq * p3;
			p3 = 66178.83658127083517939992166 + xsq * p3;
			p3 = 85145.16067533570196555001171 + xsq * p3;
			p3 = 33220.91340985722351859704442 + xsq * p3;

			double q3 = 1.0;
			q3 = 863.8367769604990967475517183 + xsq * q3;
			q3 = 37890.22974577220264142952256 + xsq * q3;
			q3 = 400294.4358226697511708610813 + xsq * q3;
			q3 = 1419460.669603720892855755253 + xsq * q3;
			q3 = 1819458.042243997298924553839 + xsq * q3;
			q3 = 708712.8194102874357377502472 + xsq * q3;

			pzero = p2 / q2;
			qzero = 8 * p3 / q3 / x;
		}
		
		/// <summary>
		/// 
		/// </summary>
		/// <param name="x"></param>
		/// <returns></returns>
		private static double Y0Hart6236(double x)
		{
			double[] P = new double[]
										{
											-.183045799017701396664E7,
											+.43903566116535873858589E7,
											-.36371585984454086979306E6,
											+ .103090170744451243873794E5,
											-.139860124593012329613322E3,
											+.103956105388106806747899E1,
											-.443191044623769135330285E-2,
											+.10431688209536418185884E-4,
											-.1082247150510209017425E-7
										};

			double[] Q = new double[]
										{									
											.2480151036050941868188E8,
											+.1970191258235426791032E6,
											+.66332162333336530231E3,
											1.0							
										};

			return PolynomialValue(P, x * x) / PolynomialValue(Q, x * x);
		}

		/// <summary>
		/// 
		/// </summary>
		/// <param name="x"></param>
		/// <param name="pzero"></param>
		/// <param name="qzero"></param>
		private static void Y0Asymptotic(double x, ref double pzero, ref double qzero)
		{
			double xsq = 64.0 / (x * x);

			double p2 = 0.0;
			p2 = 2485.271928957404011288128951 + xsq * p2;
			p2 = 153982.6532623911470917825993 + xsq * p2;
			p2 = 2016135.283049983642487182349 + xsq * p2;
			p2 = 8413041.456550439208464315611 + xsq * p2;
			p2 = 12332384.76817638145232406055 + xsq * p2;
			p2 = 5393485.083869438325262122897 + xsq * p2;

			double q2 = 1.0;
			q2 = 2615.700736920839685159081813 + xsq * q2;
			q2 = 156001.7276940030940592769933 + xsq * q2;
			q2 = 2025066.801570134013891035236 + xsq * q2;
			q2 = 8426449.050629797331554404810 + xsq * q2;
			q2 = 12338310.22786324960844856182 + xsq * q2;
			q2 = 5393485.083869438325560444960 + xsq * q2;

			double p3 = -0.0;
			p3 = -4.887199395841261531199129300 + xsq * p3;
			p3 = -226.2630641933704113967255053 + xsq * p3;
			p3 = -2365.956170779108192723612816 + xsq * p3;
			p3 = -8239.066313485606568803548860 + xsq * p3;
			p3 = -10381.41698748464093880530341 + xsq * p3;
			p3 = -3984.617357595222463506790588 + xsq * p3;
			double q3 = 1.0;
			q3 = 408.7714673983499223402830260 + xsq * q3;
			q3 = 15704.89191515395519392882766 + xsq * q3;
			q3 = 156021.3206679291652539287109 + xsq * q3;
			q3 = 533291.3634216897168722255057 + xsq * q3;
			q3 = 666745.4239319826986004038103 + xsq * q3;
			q3 = 255015.5108860942382983170882 + xsq * q3;

			pzero = p2 / q2;
			qzero = 8 * p3 / q3 / x;
		}

		
		/// <summary>
		/// 
		/// </summary>
		/// <param name="x"></param>
		/// <returns></returns>
		private static double Y1Hart6438(double x)
		{
			double[] P = new double[]
										{
											-.4900604943822124565799944536E13,
											+ .1275274390962393244694210779E13,
											-.5153438139286723656674106255E11,
											+.7349264551086296529837065763E9,
											-.4237922726264323779280112898E7,
											+.85119379358423147679302911E4					
										};

			double[] Q = new double[]
										{
											+ .2499580570086352587715671109E14,
											+.4244419664431514527800607002E12,
											+.3733650367205593523516366846E10,
											+.224590400248816548062481827E8,
											+ .1020426050426459625161062486E6,
											+.3549632885037702985804409141E3,
											1.0
										};

			return PolynomialValue(P, x * x) / PolynomialValue(Q, x * x);
		}

		/// <summary>
		/// Local value of a real polynomial by Horner.
		/// </summary>
		/// <param name="coefficients">The coefficients of the polynomial, inreasing order.</param>
		/// <param name="x">The argument.</param>
		/// <returns>The value of the polynomial</returns>
		/// <exception cref="ArgumentException">Thrown if the array is of order zero.</exception>
		private static double PolynomialValue(double[] coefficients, double x)
		{
			switch (coefficients.Length)
			{
				case 1:
					return coefficients[0];

				case 0:
					throw new ArgumentException("Polynomial is not defined");

				default:
					double d = coefficients[coefficients.Length - 1];

					for (int i = coefficients.Length - 2; i >= 0; i--)
					{
						d = d * x + coefficients[i];
					}

					return d;
			}
		}

		private static double ThetaReal(double x)
		{
			double[] b = new double[]
			{
				0.0,		0.0110486,	0.0,		-0.0000906,
				-0.0000252, -0.0000034, -0.0000006
			};

			double t = 8.0 / x;

			double s = b[6] * t + b[5];
			s = s * t + b[4];
			s = s * t + b[3];
			s = s * t + b[2];
			s = s * t + b[1];
			s = s * t + b[0];

			return s;
		}

		private static double ThetaImagine(double x)
		{
			double[] b = new double[]
			{
				-0.3926991,	-0.0110485,	-0.0009765,	-0.0000901,
				0.0,		0.0000051,	0.0000019
			};

			double t = 8.0 / x;

			double s = b[6] * t + b[5];
			s = s * t + b[4];
			s = s * t + b[3];
			s = s * t + b[2];
			s = s * t + b[1];
			s = s * t + b[0];

			return s;
		}

		/// <summary>
		/// Complexes value of I0(z) by series expansion.
		/// Equation (6.1.9) from Shanjie Zhang, Jianming Jin, Computation of special functions. ISBN 0-471-11963-6.
		/// </summary>
		/// <param name="z">The complex argument.</param>
		/// <returns>The value of I0(z).</returns>
		private static Complex ComplexI0BySeries(Complex z)
		{
			Complex result	= 0;
			Complex term	= 1;
			int n			= 0;
			Complex z2		= 0.25 * z * z;

			while (term.Magnitude > FineAccuracy)
			{
				result		+= term;
				term		*= z2;
				term		/= ((n + 1) * (n + 1));
				n++;
			}

			return result;
		}

		/// <summary>
		/// Complexes value of I0(z) by asymptotic expansion.
		/// Equation (6.2.1) from Shanjie Zhang, Jianming Jin, Computation of special functions. ISBN 0-471-11963-6.
		/// </summary>
		/// <param name="z">The complex argument.</param>
		/// <returns>The value of I0(z).</returns>
		private static Complex ComplexI0ByAsymptotic(Complex z)
		{
			// a_k = (k!!)^2 / (k! * 8^k)
			double[] coefficients	= new double[]
			{
				1,
				0.125,
				0.0703125,
				0.0732421875,
				0.112152099609375,
				0.227108001708984,
				0.572501420974731,
				1.72772750258446,
				6.07404200127348,
				24.3805296995561,
				110.017140269247,
				551.335896122021,
				3038.09051092238
			};

			int k0			= z.Magnitude >= 50.0 ? 7 : z.Magnitude >= 35.0 ? 9 : coefficients.Length - 1;

			Complex sum		= 0;
			Complex term	= 1;
			Complex zr		= 1.0 / z;

			for (int k = 0; k <= k0; k++)
			{
				sum			+= term * coefficients[k];
				term		*= zr;
			}

			return sum * Complex.Exp(z) / Complex.Sqrt(2 * Math.PI * z);
		}

		/// <summary>
		/// Complexes value of K0(z) by series expansion.
		/// Equation (6.1.11) from Shanjie Zhang, Jianming Jin, Computation of special functions. ISBN 0-471-11963-6.
		/// </summary>
		/// <param name="z">The complex argument.</param>
		/// <returns>The value of K0(z).</returns>
		private static Complex ComplexK0BySeries(Complex z)
		{
			Complex sum			= Complex.Zero;
			Complex term		= Complex.One;
			Complex argument	= z * z / 4;
			int n				= 1;
			double running		= 0.0;
			double accuracy		= Double.MaxValue;

			do
			{
				term			*= argument / (n * n);
				running			+= 1.0 / n;
				sum				+= term * running;
				n++;
				accuracy		= term.Magnitude;
			}
			while (accuracy >= MachineAccuracy);

			return -I0(z) * (Complex.Log(z/2) + Gamma) + sum;
		}

		/// <summary>
		/// Complexes value of K0(z) by asymptotic expansion.
		/// Equation (6.2.3) from Shanjie Zhang, Jianming Jin, Computation of special functions. ISBN 0-471-11963-6.
		/// </summary>
		/// <param name="z">The complex argument.</param>
		/// <returns>The value of K0(z).</returns>
		private static Complex ComplexK0ByAsymptotic(Complex z)
		{
			double[] coefficients = new double[]
			{
				0.125E0, 
				0.2109375E0, 
				1.0986328125E0, 
				1.1775970458984E01, 
				2.1461706161499E02, 
				5.9511522710323E03, 
				2.3347645606175E05, 
				7.2312234987631E07,
				8.401390346421E08,
				7.2031420482627E10
			};

			Complex z2		= z * z;
			Complex zr		= 1.0 / z2;
			Complex sum		= 1;
			Complex term	= zr;

			for (int i = 1; i < coefficients.Length; i++)
			{
				sum			+= coefficients[i] * term;
				term		*= zr;
			}

			return sum / ((2 * z) * I0(z));
		}
		#endregion

		#region Private static members		
		/// <summary>
		/// The maximum value of argument to compute I0 by series expansion
		/// (following Zhang and Jin).
		/// </summary>
		private static double	_maxSeriesArgumentI0	= 18;

		/// <summary>
		/// The maximum value of argument to compute K0 by series expansion
		/// (following Zhang and Jin).
		/// </summary>
		private static double	_maxSeriesArgumentK0	= 9;
		#endregion
	}
}
