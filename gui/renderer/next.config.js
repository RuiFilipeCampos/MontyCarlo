

  module.exports = {
    async rewrites() {
      return [
        {
          source: '/api/:path*',
          destination: 'http://127.0.0.1:5000/:path*/' // Proxy to Backend
        }
      ]
    }
  }
  