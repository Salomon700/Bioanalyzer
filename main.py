import os
from website import create_app

app = create_app()

if __name__ == '__main__':
    # Get the port from the environment variable or default to 5000
    port = int(os.environ.get('PORT', 5000))

    # Run the app on the specified host and port
    app.run(host='0.0.0.0', port=port)